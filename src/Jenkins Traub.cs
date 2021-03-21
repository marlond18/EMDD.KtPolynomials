using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using static System.Math;
using static EMDD.KtPolynomials.ComplexNumberExtensions;
using EMDD.KtPolynomials;
using KtExtensions;

namespace EMDD.KtPolynomials
{
    public static class ComplexPolynomialRootFinder
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;

        private static (Complex[] zeroRoots, Complex[] remaining) RemoveZeroRoots(Complex[] input)
        {
            var length = input.Length;
            return TestWhile(() => 0, e => e + 1, e => input[length - 1 - e] == new Complex(0, 0), e => (((Complex)0).MakeArrayFor(e), input.Take(length - e).ToArray()));
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            if (input[0] == 0) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var (zeroRoots, remaining) = RemoveZeroRoots(input);
            return zeroRoots.Concat(remaining.ScaleCoefficients().PartRoots()).ToList();
        }

        private static Complex[] ScaleCoefficients(this Complex[] input) => input.ScaleCoefficients(input.Scale());

        private static Complex[] ScaleCoefficients(this Complex[] input, double scale)
            => input.Select(coeff => coeff * scale).ToArray();

        private static IEnumerable<Complex> PartRoots(this Complex[] p)
        {
            if (p.Length == 2) yield return ComplexDivide(-p[1], p[0]);
            if (p.Length == 3) foreach (var root in Quadratic(p[0], p[1], p[2])) yield return root;
            if (p.Length <= 3) yield break;
            foreach (var root in RootsBeyond3rdDegree(p)) yield return root;
        }

        private static IEnumerable<Complex> RootsBeyond3rdDegree(Complex[] p)
        {
            var shift = new Complex(0.70710678, -0.70710678);
            var complex90 = new Complex(-0.060756474, -0.99756405);
            var scale = Cauchy(p.Select(coeff => coeff.Modulus()).ToArray());
            for (var i = 1; i <= 2; i++)
            {
                var h = NoShift(p);
                for (var j = 1; j <= 9; j++)
                {
                    shift *= complex90;
                    var (zero, converges, qp) = FxShift(10 * j, shift * scale, p, h);
                    if (!converges) continue;
                    yield return zero;
                    foreach (var root in PartRoots(qp)) yield return root;
                    yield break;
                }
            }
        }

        private static IEnumerable<Complex> Quadratic(Complex a, Complex b, Complex c)
        {
            var determ = (b * b) - (4 * a * c);
            yield return (-b + Complex.Sqrt(determ)) / (2 * a);
            yield return (-b - Complex.Sqrt(determ)) / (2 * a);
        }

        private static Complex[] Derivative(this Complex[] p)
        {
            var nn = p.Length - 1;
            var h = new Complex[nn];
            for (var i = 0; i < nn; i++) h[i] = (nn - i) * p[i] / nn;
            return h;
        }

        private static Complex[] NoShift(Complex[] p)
        {
            var nn = p.Length - 1;
            return p.Derivative().AggregateTimes(h2 => h2.Take(nn - 1).ToArray().NextH(TestKungHindiZeroComplexNumber(h2[nn - 1], p[nn - 1]), ComplexDivide(-p[nn], h2[nn - 1]), p), 5);
        }

        private static bool TestKungHindiZeroComplexNumber(Complex p1, Complex p2) => p1.Modulus() > DblEpsilon * 10 * p2.Modulus();

        private static (Complex Zero, bool Converge, Complex[] qp) FxShift(int iterations, Complex s, Complex[] p, Complex[] h)
        {
            var z = new Complex(0, 0);
            var (pv, qp) = Deflate(p, s);
            var (t, qh) = Calct(s, pv, h);
            var test = 1;
            for (var j = 1; j <= iterations; j++)
            {
                h = qh.NextH(t != new Complex(0, 0), t, qp);
                var ot = t.Clone();
                (t, qh) = Calct(s, pv, h);
                z = s + t;
                if (!(t == new Complex(0, 0) || test != 1 || j == 12) && (t - ot).Modulus() < 0.5 * z.Modulus())
                {
                    var (zero, converge, qp2) = Vrshft(z, p, h);
                    if (converge) return (zero, true, qp2);
                    test = 0;
                    (pv, qp) = Deflate(p, s);
                    (t, qh) = Calct(s, pv, h);
                }
            }
            return Vrshft(z, p, h);
        }

        private static (Complex Zero, bool converge, Complex[] qp) Vrshft(Complex z, Complex[] p, Complex[] h)
        {
            double relstp = 0;
            var b = 0;
            var s = z.Clone();
            var (pv, qp) = p.Deflate(s);
            var mp = pv.Modulus();
            if (mp <= 20 * ErrorEvaluation(qp, s.Modulus(), mp)) return (s, true, qp);
            var omp = mp;
            LastPart();
            for (int i = 2; i <= 10; i++)
            {
                (pv, qp) = p.Deflate(s);
                mp = pv.Modulus();
                if (mp <= 20 * ErrorEvaluation(qp, s.Modulus(), mp)) return (s, true, qp);
                var pumasa = b == 1 || mp < omp || relstp >= 0.05;
                if (pumasa && mp * 0.1 > omp) return (z, false, qp);
                if (!pumasa)
                {
                    b = 1;
                    s *= SMultiplier(relstp);
                    (pv, qp) = p.Deflate(s);
                    h = h.AggregateTimes(h2 => CalcTH(h2, s, qp, pv), 5);
                }
                omp = pumasa ? mp : DblMax;
                LastPart();
            }
            return (z, false, qp);
            void LastPart()
            {
                h = CalcTH(h, s, qp, pv);
                var (t1, _) = Calct(s, pv, h);
                if (t1 != new Complex(0, 0))
                {
                    relstp = t1.Modulus() / s.Modulus();
                    s += t1;
                }
            }
        }

        private static T AggregateTimes<T>(this T seed, Func<T, T> func, int times)
        {
            var sds = seed;
            for (int i = 0; i < times; i++)
            {
                sds = func(sds);
            }
            return sds;
        }

        private static Complex[] CalcTH(Complex[] h, Complex s, Complex[] qp, Complex pv)
        {
            var (t, qh) = Calct(s, pv, h);
            return qh.NextH(t != new Complex(0, 0), t, qp);
        }

        private static Complex SMultiplier(double relstp)
        {
            var r1 = Sqrt(relstp < DblEpsilon ? DblEpsilon : relstp);
            return new Complex(1 + r1, r1);
        }

        private static (Complex t, Complex[] qh) Calct(Complex s, Complex pv, Complex[] h)
        {
            var (hv, qh) = h.Deflate(s);
            return (TestKungHindiZeroComplexNumber(hv, h[^1]) ? ComplexDivide(-pv, hv) : new Complex(0, 0), qh);
        }

        private static Complex[] NextH(this Complex[] _qh, bool firstOption, Complex t, Complex[] _qp) => firstOption ? _qh.IterationNaNaman(_qp[0], (e, i) => (t * e) + _qp[i]) : _qh.IterationNaNaman(new Complex(0, 0), (e, _) => e);

        private static T[] IterationNaNaman<T>(this T[] enu, T init, Func<T, int, T> func)
        {
            var length = enu.Length;
            var h = new T[length + 1];
            h[0] = init;
            for (var j = 1; j <= length; j++)
            {
                h[j] = func(enu[j - 1], j);
            }
            return h;
        }

        private static (Complex, Complex[]) Deflate(this Complex[] p2, Complex s) => p2.Aggregate((current, prev) => (prev * s) + current).PopLast();

        private static (T last, T[] remaining) PopLast<T>(this IEnumerable<T> enu) => (enu.Last(), enu.Take(enu.Count() - 1).ToArray());

        private static double ErrorEvaluation(Complex[] q2, double ms, double mp)
        {
            var mre = 2.0 * Sqrt(2.0) * DblEpsilon;
            return (q2.Aggregate(q2[0].Modulus() * mre / (DblEpsilon + mre), (e2, elem) => (e2 * ms) + elem.Modulus()) * (DblEpsilon + mre)) - (mp * mre);
        }

        private static double Cauchy(double[] pt)
        {
            pt.NegateLast();
            return pt.GetFinalX();
        }

        private static double GetFinalX(this double[] pt2) => TestWhile(() => pt2.GetXdX(pt2.GetSecondX()), d => pt2.GetXdX(d.x), d => d.ratio > 0.005, d => d.x);

        private static T TestWhile<T, T2>(Func<T2> start, Func<T2, T2> body, Func<T2, bool> test, Func<T2, T> @return)
        {
            var init = start();
            while (test(init))
            {
                init = body(init);
            }
            return @return(init);
        }

        private static (double x, double ratio) GetXdX(this double[] pt2, double x)
        {
            var dx = Getdx(pt2, x);
            var x2 = x - dx;
            return (x2, Abs(dx / x2));
        }

        private static double Getdx(double[] pt2, double x)
        {
            var (f1, q) = pt2.Deflate(x);
            var (f2, _) = q.Deflate(x);
            return f1 / f2;
        }

        private static void NegateLast(this double[] pt2)
        {
            var nn = pt2.Length - 1;
            pt2[nn] = -pt2[nn];
        }

        private static double GetSecondX(this double[] pt2) => TestWhile(() => GetInitialX(pt2) * 0.1, d => d * 0.1, d => pt2.Deflate(d).result > 0, d => d);

        private static (double result, double[] remaining) Deflate(this double[] coeff, double s) => coeff.Aggregate((current, prev) => (prev * s) + current).PopLast();

        private static T[] Aggregate<T>(this T[] enu, Func<T, T, T> func)
        {
            if (func == null) return null;
            var length = enu.Length;
            var q = enu.ToArray();
            for (int i = 1; i < length; i++)
            {
                q[i] = func(enu[i], q[i - 1]);
            }
            return q;
        }

        private static double GetInitialX(double[] pt2)
        {
            var nn = pt2.Length - 1;
            var x = Exp(Log(-pt2[nn] / pt2[0])) / nn;
            if (!pt2[nn - 1].NearZero())
            {
                var xm = -pt2[nn] / pt2[nn - 1];
                if (xm < x) return xm;
            }
            return x;
        }

        private static double Scale(this Complex[] pt)
        {
            const double lo = DblMin / DblEpsilon;
            var (min, max) = pt.GetMinMaxModulus();
            if (min >= lo && max <= Sqrt(DblMax)) return 1;
            var x = lo / min;
            return (2.0).RaiseTo(Round(Log(x <= 1 ? 1 / (Sqrt(max) * Sqrt(min)) : DblMax / x > max ? 1 : x, 2) + 0.5));
        }

        private static (double min, double max) GetMinMaxModulus(this Complex[] pt)
        {
            var max = 0.0;
            var min = DblMax;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                var x = pt[i].Modulus();
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            return (min, max);
        }
    }

    internal static class ComplexNumberExtensions
    {
        private const double DblMax = 1.79769313486232E+307;

#pragma warning disable RCS1224 // Make method an extension method.
        public static Complex ComplexDivide(Complex a, Complex b) =>
#pragma warning restore RCS1224 // Make method an extension method.
            b.Real.NearZero() && b.Imaginary.NearZero() ?
            new Complex(DblMax, DblMax) : a / b;

        public static double Modulus(this Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai) return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            if (ar > ai) return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            return ar * Sqrt(2.0);
        }

        internal static Complex Clone(this Complex val) => new(val.Real, val.Imaginary);
    }
}