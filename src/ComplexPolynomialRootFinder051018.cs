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
    public static class ComplexPolynomialRootFinder051018
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

        private static Complex[] ScaleCoefficients(this Complex[] input)
        {
            var scale = input.Scale();
            return input.Select(coeff => coeff * scale).ToArray();
        }

        private static IEnumerable<Complex> PartRoots(this Complex[] p)
        {
            if (p.Length == 2)
            {
                yield return ComplexDivide(-p[1], p[0]);
            }
            if (p.Length == 3)
            {
                foreach (var root in Quadratic(p[0], p[1], p[2]))
                {
                    yield return root;
                }
            }
            if (p.Length <= 3) yield break;
            var shift = new Complex(0.70710678, -0.70710678);
            var complex90 = new Complex(-0.060756474, -0.99756405);
            var scale = Cauchy(p.Select(ComplexNumberExtensions.Modulus).ToArray());
            for (var i = 1; i <= 2; i++)
            {
                var h = NoShift(p);
                for (var j = 1; j <= 9; j++)
                {
                    shift *= complex90;
                    var (zero, converges, qp) = FxShift(10 * j, shift * scale, p, h);
                    if (converges)
                    {
                        yield return zero;
                        foreach (var root in PartRoots(qp))
                        {
                            yield return root;
                        }
                        yield break;
                    }
                }
            }
        }

        private static IEnumerable<Complex> Quadratic(Complex a, Complex b, Complex c)
        {
            var determ = (b * b) - (4 * a * c);
            yield return (-b + Complex.Sqrt(determ)) / (2 * a);
            yield return (-b - Complex.Sqrt(determ)) / (2 * a);
        }

        private static Complex[] Derivative(Complex[] p)
        {
            var nn = p.Length - 1;
            var h = new Complex[nn];
            for (var i = 0; i < nn; i++)
            {
                h[i] = (nn - i) * p[i] / nn;
            }
            return h;
        }

        private static Complex[] NoShift(Complex[] p)
        {
            var nn = p.Length - 1;
            var h = Derivative(p);
            for (var jj = 1; jj <= 5; jj++)
            {
                h = h.Take(nn - 1).ToArray().NextH(h[nn - 1].Modulus() > DblEpsilon * 10 * p[nn - 1].Modulus(), ComplexDivide(-p[nn], h[nn - 1]), p);
                //    if (ComplexModulus(h[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                //    {
                //        var t = ComplexDivide(-p[nn], h[nn - 1]);
                //        for (int i2 = h.Length - 1; i2 >= 1; i2--)
                //        {
                //            h[i2] = (t * h[i2 - 1]) + p[i2];
                //        }
                //        h[0] = p[0].Clone();
                //    }
                //    else
                //    {
                //        for (int i2 = h.Length - 1; i2 >= 1; i2--)
                //        {
                //            h[i2] = h[i2 - 1];
                //        }
                //        h[0] = new Complex(0, 0);
                //    }
            }
            return h;
        }

        private static (Complex Zero, bool Converge, Complex[] qp) FxShift(int iterations, Complex s, Complex[] p, Complex[] h)
        {
            var z = new Complex(0, 0);
            var (pv, qp) = Deflate(p, s);
            var (t, qh) = Calct(s, pv, h);
            var test = 1;
            var pasd = 0;
            for (var j = 1; j <= iterations; j++)
            {
                h = qh.NextH(t != new Complex(0, 0), t, qp);
                var ot = t.Clone();
                (t, qh) = Calct(s, pv, h);
                z = s + t;
                var pumasa = !(t == new Complex(0, 0) || test != 1 || j == 12) && (t - ot).Modulus() < 0.5 * z.Modulus();
                if (pumasa && pasd == 1)
                {
                    var (zero, converge, qp2) = Vrshft(z, p, h);
                    if (converge) return (zero, true, qp2);
                    test = 0;
                    (pv, qp) = Deflate(p, s);
                    (t, qh) = Calct(s, pv, h);
                    continue;
                }
                pasd = pumasa ? 1 : 0;
            }
            return Vrshft(z, p, h);
        }

        private static (Complex Zero, bool converge, Complex[] qp) Vrshft(Complex z, Complex[] p, Complex[] h)
        {
            double relstp = 0;
            var b = 0;
            var s = z.Clone();
            _ = new Complex[p.Length];
            Complex pv;
            Complex[] qp;
            (pv, qp) = p.Deflate(s);
            var mp = pv.Modulus();
            if (mp <= 20 * ErrorEvaluation(qp, s.Modulus(), mp)) return (s, true, qp);
            double omp = mp;
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
                    for (var j = 1; j <= 5; j++)
                    {
                        h = CalcTH(h, s, qp, pv);
                    }
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
            return (!(hv.Modulus() <= DblEpsilon * 10 * h[^1].Modulus()) ? ComplexDivide(-pv, hv) : new Complex(0, 0), qh);
        }

        private static Complex[] NextH(this Complex[] _qh, bool firstOption, Complex t, Complex[] _qp) => firstOption /*(t != new Complex(0, 0))*/ ? _qh.Tatatatata(_qp[0], (e, i) => (t * e) + _qp[i]) : _qh.Tatatatata(new Complex(0, 0), (e, _) => e);

        private static T[] Tatatatata<T>(this T[] enu, T init, Func<T, int, T> func)
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
            var (min, max) = pt.GetMinMax();
            if (min >= lo && max <= Sqrt(DblMax)) return 1;
            var x = lo / min;
            return (2.0).RaiseTo(Round(Log(x <= 1 ? 1 / (Sqrt(max) * Sqrt(min)) : DblMax / x > max ? 1 : x, 2) + 0.5));
        }

        private static (double min, double max) GetMinMax(this Complex[] pt)
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
}