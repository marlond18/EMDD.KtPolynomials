using KtExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using static System.Math;
using static System.Convert;

namespace EMDD.KtPolynomials
{
    public static class ComplexPolynomialRootFinder1ChangeNoSwift
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            var result = new List<Complex>();
            var zeroRoots = CountNumbersZeroRoots(input);
            var copyOfCoeff = new Complex[input.Length - zeroRoots];
            for (var i = 0; i < zeroRoots; i++) result.Add(new Complex(0, 0));
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var degree = copyOfCoeff.Length - 1;
            var p = new Complex[degree + 1];
            var h = new Complex[degree + 1];
            var _qp = new Complex[degree + 1];
            var _qh = new Complex[degree + 1];
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (copyOfCoeff[0] == new Complex(0, 0)) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var scale = Scale(copyOfCoeff);
            var shift = new Complex(0.70710678, -0.70710678);
            p = copyOfCoeff.Select(coeff => coeff * scale).ToArray();
            var nn = degree;
            search:
            if (nn <= 1)
            {
                result.Add(ComplexDivide(-p[1], p[0]));
                return result;
            }
            scale = Cauchy(nn, p.Select(ComplexModulus).ToArray());
            for (var cnt1 = 1; cnt1 <= 2; cnt1++)
            {
                h = NoShift(5, p);
                for (var cnt2 = 1; cnt2 <= 9; cnt2++)
                {
                    shift *= complex90;
                    var (Zero, Converge) = FxShift(nn, 10 * cnt2, shift * scale, p, h, _qp, _qh);
                    if (Converge)
                    {
                        result.Add(Zero.Clone());
                        nn--;
                        p = new Complex[nn + 1];
                        for (var i = 0; i <= nn; i++)
                        {
                            p[i] = _qp[i].Clone();
                        }
                        goto search;
                    }
                }
            }
            return result;
        }

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var tempH = new Complex[nn + 1];
            for (var i = 0; i <= nn; i++)
            {
                double xni = nn - i;
                tempH[i] = xni * p[i] / nn;
            }
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(tempH[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], tempH[nn - 1]);
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = (t * tempH[nn - i - 2]) + p[nn - i - 1];
                    }
                    tempH[0] = p[0].Clone();
                }
                else
                {
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = tempH[nn - i - 2];
                    }
                    tempH[0] = new Complex(0, 0);
                }
            }
            return tempH;
        }

        // COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
        // INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
        // APPROXIMATE ZERO IF SUCCESSFUL.
        // L2 - LIMIT OF FIXED SHIFT STEPS
        // ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
        // CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
        private static (Complex Zero, bool Converge) FxShift(int nn, int l2, Complex s, Complex[] p, Complex[] initH, Complex[] _qp, Complex[] _qh)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var sh = new Complex[p.Length];
            var z2 = new Complex(0, 0);
            var n = nn;
            var pv = PolyEvaluation(nn, s, p, _qp);
            var test = 1;
            var pasd = 0;
            // Calculate first T = -P(S)/H(S)
            var calcT = Calct(nn, s, pv, tempH, _qh);
            // Main loop for second stage
            (Complex Zero, bool converge, Complex[] H) vrshift;
            for (var j = 1; j <= l2; j++)
            {
                var ot = calcT.T.Clone();
                // Compute the next H Polynomial and new t
                tempH = NextH(nn, calcT.HisZero, calcT.T, tempH, _qh, _qp);
                calcT = Calct(nn, s, pv, tempH, _qh);
                z2 = s + calcT.T;
                // Test for convergence unless stage 3 has failed once or this
                // is the last H Polynomial
                if (!(calcT.HisZero || test != 1 || j == 12))
                {
                    if (ComplexModulus(calcT.T - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            // The weak convergence test has been passwed twice, start the third stage
                            // Iteration, after saving the current H polynomial and shift
                            for (var i = 0; i <= n - 1; i++)
                            {
                                sh[i] = tempH[i].Clone();
                            }
                            vrshift = Vrshft(nn, 10, z2, p, tempH, _qp, _qh);
                            tempH = vrshift.H;
                            z2 = vrshift.Zero;
                            if (vrshift.converge) return (z2, true);
                            //The iteration failed to converge. Turn off testing and restore h,s,pv and T
                            test = 0;
                            for (var i = 0; i <= n - 1; i++)
                            {
                                tempH[i] = sh[i].Clone();
                            }
                            pv = PolyEvaluation(nn, s, p, _qp);
                            calcT = Calct(nn, s, pv, tempH, _qh);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            // Attempt an iteration with final H polynomial from second stage
            vrshift = Vrshft(nn, 10, z2, p, tempH, _qp, _qh);
            z2 = vrshift.Zero;
            return (z2, vrshift.converge);
        }

        // CARRIES OUT THE THIRD STAGE ITERATION.
        // L3 - LIMIT OF STEPS IN STAGE 3.
        // ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
        //           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
        // CONV    -  .TRUE. IF ITERATION CONVERGES
        //
        private static (Complex Zero, bool converge, Complex[] H) Vrshft(int nn, int l3, Complex z2, Complex[] p, Complex[] initH, Complex[] _qp, Complex[] _qh)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            // Main loop for stage three
            for (int i = 1; i <= l3; i++)
            {
                // Evaluate P at S and test for convergence
                var pv = PolyEvaluation(nn, s, p, _qp);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(nn, _qp, ms, mp, DblEpsilon, 2.0 * Sqrt(2.0) * DblEpsilon))
                {
                    z2 = s;
                    return (z2, true, tempH);
                }
                (bool HisZero, Complex T) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed
                        // shift steps into the cluster to force one zero to dominate
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        pv = PolyEvaluation(nn, s, p, _qp);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(nn, s, pv, tempH, _qh);
                            tempH = NextH(nn, calcT.HisZero, calcT.T, tempH, _qh, _qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }

                    // Exit if polynomial value increase significantly
                    if (mp * 0.1 > omp) return (z2, false, tempH);
                }
                omp = mp;
                _20:
                // Calculate next iterate
                calcT = Calct(nn, s, pv, tempH, _qh);
                tempH = NextH(nn, calcT.HisZero, calcT.T, tempH, _qh, _qp);
                calcT = Calct(nn, s, pv, tempH, _qh);
                if (!calcT.HisZero)
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false, tempH);
        }

        // COMPUTES  T = -P(S)/H(S).
        // BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
        private static (bool HisZero, Complex T) Calct(int nn, Complex s, Complex pv, Complex[] initH, Complex[] _qh)
        {
            var hv = PolyEvaluation(nn - 1, s, initH, _qh);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn - 1]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (bol, t);
        }

        // CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
        // BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
        private static Complex[] NextH(int nn, bool bol, Complex t, Complex[] initH, Complex[] _qh, Complex[] _qp)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            if (!bol)
            {
                for (var j = 1; j <= nn - 1; j++)
                {
                    tempH[j] = (t * _qh[j - 1]) + _qp[j];
                }
                tempH[0] = _qp[0];
                return tempH;
            }
            for (var j = 1; j <= nn - 1; j++)
            {
                tempH[j] = _qh[j - 1];
            }
            tempH[0] = new Complex(0, 0);
            return tempH;
        }
        // EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
        // PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
        private static Complex PolyEvaluation(int nn2, Complex s2, Complex[] p2, Complex[] q2)
        {
            q2[0] = p2[0].Clone();
            var pv = q2[0].Clone();
            for (var i = 1; i <= nn2; i++)
            {
                pv = (pv * s2) + p2[i];
                q2[i] = pv.Clone();
            }
            return pv;
        }

        private static double ErrorEvaluation(int nn2, Complex[] q2, double ms, double mp, double are, double mre)
        {
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= nn2; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        // CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
        // POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
        //
        private static double Cauchy(int nn2, double[] pt2)
        {
            double xm;
            double f;
            pt2[nn2] = -pt2[nn2];
            // Compute upper estimate bound
            var n = nn2;
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / n;
            if (!pt2[n - 1].NearZero())
            {
                xm = -pt2[nn2] / pt2[n - 1];
                if (xm < x) x = xm;
            }
            // Chop the interval (0,x) until f < 0
            while (true)
            {
                xm = x * 0.1;
                f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            // Do Newton iteration until x converges to two decimal places
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= n - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            // Find largest and smallest moduli of coefficients
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            // Scale only if there are very large or very small components
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        // COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
        private static Complex ComplexDivide(Complex a, Complex b)
        {
            if (b.Real.NearZero() && b.Imaginary.NearZero())
            {
                return new Complex(DblMax, DblMax);
            }
            return a / b;
        }

        // MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }

    public static class ComplexPolynomialRootFinder2PartialNNRemoval
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            var result = new List<Complex>();
            var zeroRoots = CountNumbersZeroRoots(input);
            var copyOfCoeff = new Complex[input.Length - zeroRoots];
            for (var i = 0; i < zeroRoots; i++) result.Add(new Complex(0, 0));
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var degree = copyOfCoeff.Length - 1;
            var _qp = new Complex[degree + 1];
            var _qh = new Complex[degree + 1];
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (copyOfCoeff[0] == new Complex(0, 0)) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var scale = Scale(copyOfCoeff);
            var shift = new Complex(0.70710678, -0.70710678);
            var p = copyOfCoeff.Select(coeff => coeff * scale).ToArray();
            var nn = degree;
            search:
            if (nn <= 1)
            {
                result.Add(ComplexDivide(-p[1], p[0]));
                return result;
            }
            scale = Cauchy(p.Select(ComplexModulus).ToArray());
            for (var cnt1 = 1; cnt1 <= 2; cnt1++)
            {
                var h = NoShift(5, p);
                for (var cnt2 = 1; cnt2 <= 9; cnt2++)
                {
                    shift *= complex90;
                    var (Zero, Converge) = FxShift(10 * cnt2, shift * scale, p, h, _qp, _qh);
                    if (Converge)
                    {
                        result.Add(Zero.Clone());
                        nn--;
                        p = new Complex[nn + 1];
                        for (var i = 0; i <= nn; i++)
                        {
                            p[i] = _qp[i].Clone();
                        }
                        goto search;
                    }
                }
            }
            return result;
        }

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var tempH = new Complex[nn];
            for (var i = 0; i < nn; i++)
            {
                double xni = nn - i;
                tempH[i] = xni * p[i] / nn;
            }
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(tempH[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], tempH[nn - 1]);
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = (t * tempH[nn - i - 2]) + p[nn - i - 1];
                    }
                    tempH[0] = p[0].Clone();
                }
                else
                {
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = tempH[nn - i - 2];
                    }
                    tempH[0] = new Complex(0, 0);
                }
            }
            return tempH;
        }

        // COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
        // INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
        // APPROXIMATE ZERO IF SUCCESSFUL.
        // L2 - LIMIT OF FIXED SHIFT STEPS
        // ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
        // CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
        private static (Complex Zero, bool Converge) FxShift(int l2, Complex s, Complex[] p, Complex[] initH, Complex[] _qp, Complex[] _qh)
        {
            var nn = p.Length - 1;
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var sh = new Complex[p.Length];
            var z2 = new Complex(0, 0);
            var n = nn;
            var pv = PolyEvaluation(nn, s, p, _qp);
            var test = 1;
            var pasd = 0;
            // Calculate first T = -P(S)/H(S)
            var calcT = Calct(s, pv, tempH, _qh);
            // Main loop for second stage
            (Complex Zero, bool converge, Complex[] H) vrshift;
            for (var j = 1; j <= l2; j++)
            {
                var ot = calcT.T.Clone();
                // Compute the next H Polynomial and new t
                tempH = NextH(calcT.HisZero, calcT.T, tempH, _qh, _qp);
                calcT = Calct(s, pv, tempH, _qh);
                z2 = s + calcT.T;
                // Test for convergence unless stage 3 has failed once or this
                // is the last H Polynomial
                if (!(calcT.HisZero || test != 1 || j == 12))
                {
                    if (ComplexModulus(calcT.T - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            // The weak convergence test has been passwed twice, start the third stage
                            // Iteration, after saving the current H polynomial and shift
                            for (var i = 0; i <= n - 1; i++)
                            {
                                sh[i] = tempH[i].Clone();
                            }
                            vrshift = Vrshft(10, z2, p, tempH, _qp, _qh);
                            tempH = vrshift.H;
                            z2 = vrshift.Zero;
                            if (vrshift.converge) return (z2, true);
                            //The iteration failed to converge. Turn off testing and restore h,s,pv and T
                            test = 0;
                            for (var i = 0; i <= n - 1; i++)
                            {
                                tempH[i] = sh[i].Clone();
                            }
                            pv = PolyEvaluation(nn, s, p, _qp);
                            calcT = Calct(s, pv, tempH, _qh);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            // Attempt an iteration with final H polynomial from second stage
            vrshift = Vrshft(10, z2, p, tempH, _qp, _qh);
            z2 = vrshift.Zero;
            return (z2, vrshift.converge);
        }

        // CARRIES OUT THE THIRD STAGE ITERATION.
        // L3 - LIMIT OF STEPS IN STAGE 3.
        // ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
        //           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
        // CONV    -  .TRUE. IF ITERATION CONVERGES
        //
        private static (Complex Zero, bool converge, Complex[] H) Vrshft(int l3, Complex z2, Complex[] p, Complex[] initH, Complex[] _qp, Complex[] _qh)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var nn = p.Length - 1;
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            // Main loop for stage three
            for (int i = 1; i <= l3; i++)
            {
                // Evaluate P at S and test for convergence
                var pv = PolyEvaluation(nn, s, p, _qp);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(nn, _qp, ms, mp, DblEpsilon, 2.0 * Sqrt(2.0) * DblEpsilon))
                {
                    z2 = s;
                    return (z2, true, tempH);
                }
                (bool HisZero, Complex T) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed
                        // shift steps into the cluster to force one zero to dominate
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        pv = PolyEvaluation(nn, s, p, _qp);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(s, pv, tempH, _qh);
                            tempH = NextH(calcT.HisZero, calcT.T, tempH, _qh, _qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }

                    // Exit if polynomial value increase significantly
                    if (mp * 0.1 > omp) return (z2, false, tempH);
                }
                omp = mp;
                _20:
                // Calculate next iterate
                calcT = Calct(s, pv, tempH, _qh);
                tempH = NextH(calcT.HisZero, calcT.T, tempH, _qh, _qp);
                calcT = Calct(s, pv, tempH, _qh);
                if (!calcT.HisZero)
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false, tempH);
        }

        // COMPUTES  T = -P(S)/H(S).
        // BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
        private static (bool HisZero, Complex T) Calct(Complex s, Complex pv, Complex[] initH, Complex[] _qh)
        {
            var nn = initH.Length - 1;
            var hv = PolyEvaluation(nn, s, initH, _qh);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (bol, t);
        }

        // CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
        // BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
        private static Complex[] NextH(bool bol, Complex t, Complex[] initH, Complex[] _qh, Complex[] _qp)
        {
            var nn = initH.Length;
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            if (!bol)
            {
                for (var j = 1; j <= nn - 1; j++)
                {
                    tempH[j] = (t * _qh[j - 1]) + _qp[j];
                }
                tempH[0] = _qp[0];
                return tempH;
            }
            for (var j = 1; j <= nn - 1; j++)
            {
                tempH[j] = _qh[j - 1];
            }
            tempH[0] = new Complex(0, 0);
            return tempH;
        }
        // EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
        // PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
        private static Complex PolyEvaluation(int nn2, Complex s2, Complex[] p2, Complex[] q2)
        {
            q2[0] = p2[0].Clone();
            var pv = q2[0].Clone();
            for (var i = 1; i <= nn2; i++)
            {
                pv = (pv * s2) + p2[i];
                q2[i] = pv.Clone();
            }
            return pv;
        }

        private static double ErrorEvaluation(int nn2, Complex[] q2, double ms, double mp, double are, double mre)
        {
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= nn2; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        // CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
        // POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
        //
        private static double Cauchy(double[] pt2)
        {
            double xm;
            double f;
            var nn2 = pt2.Length - 1;
            pt2[nn2] = -pt2[nn2];
            // Compute upper estimate bound
            var n = nn2;
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / n;
            if (!pt2[n - 1].NearZero())
            {
                xm = -pt2[nn2] / pt2[n - 1];
                if (xm < x) x = xm;
            }
            // Chop the interval (0,x) until f < 0
            while (true)
            {
                xm = x * 0.1;
                f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            // Do Newton iteration until x converges to two decimal places
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= n - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            // Find largest and smallest moduli of coefficients
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            // Scale only if there are very large or very small components
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        // COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
        private static Complex ComplexDivide(Complex a, Complex b)
        {
            if (b.Real.NearZero() && b.Imaginary.NearZero())
            {
                return new Complex(DblMax, DblMax);
            }
            return a / b;
        }

        // MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }

    public static class ComplexPolynomialRootFinder3
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            var result = new List<Complex>();
            var zeroRoots = CountNumbersZeroRoots(input);
            var copyOfCoeff = new Complex[input.Length - zeroRoots];
            for (var i = 0; i < zeroRoots; i++) result.Add(new Complex(0, 0));
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var degree = copyOfCoeff.Length - 1;
            var _qp = new Complex[degree + 1];
            var _qh = new Complex[degree + 1];
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (copyOfCoeff[0] == new Complex(0, 0)) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var scale = Scale(copyOfCoeff);
            var shift = new Complex(0.70710678, -0.70710678);
            var p = copyOfCoeff.Select(coeff => coeff * scale).ToArray();
            var nn = degree;
            search:
            if (nn <= 1)
            {
                result.Add(ComplexDivide(-p[1], p[0]));
                return result;
            }
            scale = Cauchy(p.Select(ComplexModulus).ToArray());
            for (var cnt1 = 1; cnt1 <= 2; cnt1++)
            {
                var h = NoShift(5, p);
                for (var cnt2 = 1; cnt2 <= 9; cnt2++)
                {
                    shift *= complex90;
                    var (Zero, Converge) = FxShift(10 * cnt2, shift * scale, p, h, _qp, _qh);
                    if (Converge)
                    {
                        result.Add(Zero.Clone());
                        nn--;
                        p = new Complex[nn + 1];
                        for (var i = 0; i <= nn; i++)
                        {
                            p[i] = _qp[i].Clone();
                        }
                        goto search;
                    }
                }
            }
            return result;
        }

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var tempH = new Complex[nn];
            for (var i = 0; i < nn; i++)
            {
                double xni = nn - i;
                tempH[i] = xni * p[i] / nn;
            }
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(tempH[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], tempH[nn - 1]);
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = (t * tempH[nn - i - 2]) + p[nn - i - 1];
                    }
                    tempH[0] = p[0].Clone();
                }
                else
                {
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = tempH[nn - i - 2];
                    }
                    tempH[0] = new Complex(0, 0);
                }
            }
            return tempH;
        }

        private static (Complex Zero, bool Converge) FxShift(int l2, Complex s, Complex[] p, Complex[] initH, Complex[] _qp, Complex[] _qh)
        {
            var nn = p.Length - 1;
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var sh = new Complex[p.Length];
            var z2 = new Complex(0, 0);
            var n = nn;
            var pv = PolyEvaluation(s, p, _qp);
            var test = 1;
            var pasd = 0;
            var calcT = Calct(s, pv, tempH, _qh);
            (Complex Zero, bool converge) vrshift;
            for (var j = 1; j <= l2; j++)
            {
                var ot = calcT.T.Clone();
                tempH = NextH(calcT.HisZero, calcT.T, tempH, _qh, _qp);
                calcT = Calct(s, pv, tempH, _qh);
                z2 = s + calcT.T;
                if (!(calcT.HisZero || test != 1 || j == 12))
                {
                    if (ComplexModulus(calcT.T - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            vrshift = Vrshft(10, z2, p, tempH, _qp, _qh);
                            z2 = vrshift.Zero;
                            if (vrshift.converge) return (z2, true);
                            test = 0;
                            pv = PolyEvaluation(s, p, _qp);
                            calcT = Calct(s, pv, tempH, _qh);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            vrshift = Vrshft(10, z2, p, tempH, _qp, _qh);
            z2 = vrshift.Zero;
            return (z2, vrshift.converge);
        }

        private static (Complex Zero, bool converge) Vrshft(int l3, Complex z2, Complex[] p, Complex[] initH, Complex[] _qp, Complex[] _qh)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var nn = p.Length - 1;
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            for (int i = 1; i <= l3; i++)
            {
                var pv = PolyEvaluation(s, p, _qp);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(nn, _qp, ms, mp, DblEpsilon, 2.0 * Sqrt(2.0) * DblEpsilon))
                {
                    z2 = s;
                    return (z2, true);
                }
                (bool HisZero, Complex T) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        pv = PolyEvaluation(s, p, _qp);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(s, pv, tempH, _qh);
                            tempH = NextH(calcT.HisZero, calcT.T, tempH, _qh, _qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }
                    if (mp * 0.1 > omp) return (z2, false);
                }
                omp = mp;
                _20:
                calcT = Calct(s, pv, tempH, _qh);
                tempH = NextH(calcT.HisZero, calcT.T, tempH, _qh, _qp);
                calcT = Calct(s, pv, tempH, _qh);
                if (!calcT.HisZero)
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false);
        }

        private static (bool HisZero, Complex T) Calct(Complex s, Complex pv, Complex[] initH, Complex[] _qh)
        {
            var nn = initH.Length - 1;
            var hv = PolyEvaluation(s, initH, _qh);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (bol, t);
        }

        private static Complex[] NextH(bool bol, Complex t, Complex[] initH, Complex[] _qh, Complex[] _qp)
        {
            var nn = initH.Length;
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            if (!bol)
            {
                for (var j = 1; j <= nn - 1; j++)
                {
                    tempH[j] = (t * _qh[j - 1]) + _qp[j];
                }
                tempH[0] = _qp[0];
                return tempH;
            }
            for (var j = 1; j <= nn - 1; j++)
            {
                tempH[j] = _qh[j - 1];
            }
            tempH[0] = new Complex(0, 0);
            return tempH;
        }

        private static Complex PolyEvaluation(Complex s2, Complex[] p2, Complex[] q2)
        {
            var nn2 = p2.Length - 1;
            q2[0] = p2[0].Clone();
            var pv = q2[0].Clone();
            for (var i = 1; i <= nn2; i++)
            {
                pv = (pv * s2) + p2[i];
                q2[i] = pv.Clone();
            }
            return pv;
        }

        private static double ErrorEvaluation(int nn2, Complex[] q2, double ms, double mp, double are, double mre)
        {
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= nn2; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        private static double Cauchy(double[] pt2)
        {
            double xm;
            double f;
            var nn2 = pt2.Length - 1;
            pt2[nn2] = -pt2[nn2];
            var n = nn2;
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / n;
            if (!pt2[n - 1].NearZero())
            {
                xm = -pt2[nn2] / pt2[n - 1];
                if (xm < x) x = xm;
            }
            while (true)
            {
                xm = x * 0.1;
                f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= n - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        private static Complex ComplexDivide(Complex a, Complex b)
        {
            if (b.Real.NearZero() && b.Imaginary.NearZero())
            {
                return new Complex(DblMax, DblMax);
            }
            return a / b;
        }

        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }

    public static class ComplexPolynomialRootFinder4
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            var result = new List<Complex>();
            var zeroRoots = CountNumbersZeroRoots(input);
            var copyOfCoeff = new Complex[input.Length - zeroRoots];
            for (var i = 0; i < zeroRoots; i++) result.Add(new Complex(0, 0));
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var degree = copyOfCoeff.Length - 1;
            var _qp = new Complex[degree + 1];
            var _qh = new Complex[degree + 1];
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (copyOfCoeff[0] == new Complex(0, 0)) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var scale = Scale(copyOfCoeff);
            var shift = new Complex(0.70710678, -0.70710678);
            var p = copyOfCoeff.Select(coeff => coeff * scale).ToArray();
            var nn = degree;
            search:
            if (nn <= 1)
            {
                result.Add(ComplexDivide(-p[1], p[0]));
                return result;
            }
            scale = Cauchy(p.Select(ComplexModulus).ToArray());
            for (var cnt1 = 1; cnt1 <= 2; cnt1++)
            {
                var h = NoShift(5, p);
                for (var cnt2 = 1; cnt2 <= 9; cnt2++)
                {
                    shift *= complex90;
                    var (Zero, Converge) = FxShift(10 * cnt2, shift * scale, p, h, _qp);
                    if (Converge)
                    {
                        result.Add(Zero.Clone());
                        nn--;
                        p = new Complex[nn + 1];
                        for (var i = 0; i <= nn; i++)
                        {
                            p[i] = _qp[i].Clone();
                        }
                        goto search;
                    }
                }
            }
            return result;
        }

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var tempH = new Complex[nn];
            for (var i = 0; i < nn; i++)
            {
                double xni = nn - i;
                tempH[i] = xni * p[i] / nn;
            }
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(tempH[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], tempH[nn - 1]);
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = (t * tempH[nn - i - 2]) + p[nn - i - 1];
                    }
                    tempH[0] = p[0].Clone();
                }
                else
                {
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        tempH[nn - i - 1] = tempH[nn - i - 2];
                    }
                    tempH[0] = new Complex(0, 0);
                }
            }
            return tempH;
        }

        private static (Complex Zero, bool Converge) FxShift(int l2, Complex s, Complex[] p, Complex[] initH, Complex[] _qp)
        {
            var nn = p.Length - 1;
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var sh = new Complex[p.Length];
            var z2 = new Complex(0, 0);
            var n = nn;
            var pv = PolyEvaluation(s, p, _qp);
            var test = 1;
            var pasd = 0;
            var calcT = Calct(s, pv, tempH);
            (Complex Zero, bool converge) vrshift;
            for (var j = 1; j <= l2; j++)
            {
                var ot = calcT.T.Clone();
                tempH = NextH(calcT.HisZero, calcT.T, tempH, calcT._qh, _qp);
                calcT = Calct(s, pv, tempH);
                z2 = s + calcT.T;
                if (!(calcT.HisZero || test != 1 || j == 12))
                {
                    if (ComplexModulus(calcT.T - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            vrshift = Vrshft(10, z2, p, tempH, _qp);
                            z2 = vrshift.Zero;
                            if (vrshift.converge) return (z2, true);
                            test = 0;
                            pv = PolyEvaluation(s, p, _qp);
                            calcT = Calct(s, pv, tempH);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            vrshift = Vrshft(10, z2, p, tempH, _qp);
            z2 = vrshift.Zero;
            return (z2, vrshift.converge);
        }

        private static (Complex Zero, bool converge) Vrshft(int l3, Complex z2, Complex[] p, Complex[] initH, Complex[] _qp)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var nn = p.Length - 1;
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            for (int i = 1; i <= l3; i++)
            {
                var pv = PolyEvaluation(s, p, _qp);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(nn, _qp, ms, mp, DblEpsilon, 2.0 * Sqrt(2.0) * DblEpsilon))
                {
                    z2 = s;
                    return (z2, true);
                }
                (bool HisZero, Complex T, Complex[] qh) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        pv = PolyEvaluation(s, p, _qp);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(s, pv, tempH);
                            tempH = NextH(calcT.HisZero, calcT.T, tempH, calcT.qh, _qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }
                    if (mp * 0.1 > omp) return (z2, false);
                }
                omp = mp;
                _20:
                calcT = Calct(s, pv, tempH);
                tempH = NextH(calcT.HisZero, calcT.T, tempH, calcT.qh, _qp);
                calcT = Calct(s, pv, tempH);
                if (!calcT.HisZero)
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false);
        }

        private static (bool HisZero, Complex T, Complex[] _qh) Calct(Complex s, Complex pv, Complex[] initH)
        {
            var nn = initH.Length - 1;
            var _qh = new Complex[nn + 1];
            var hv = PolyEvaluation(s, initH, _qh);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (bol, t, _qh);
        }

        private static Complex[] NextH(bool bol, Complex t, Complex[] initH, Complex[] _qh, Complex[] _qp)
        {
            var nn = initH.Length;
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            if (!bol)
            {
                for (var j = 1; j <= nn - 1; j++)
                {
                    tempH[j] = (t * _qh[j - 1]) + _qp[j];
                }
                tempH[0] = _qp[0];
                return tempH;
            }
            for (var j = 1; j <= nn - 1; j++)
            {
                tempH[j] = _qh[j - 1];
            }
            tempH[0] = new Complex(0, 0);
            return tempH;
        }

        private static Complex PolyEvaluation(Complex s2, Complex[] p2, Complex[] q2)
        {
            var nn2 = p2.Length - 1;
            q2[0] = p2[0].Clone();
            var pv = q2[0].Clone();
            for (var i = 1; i <= nn2; i++)
            {
                pv = (pv * s2) + p2[i];
                q2[i] = pv.Clone();
            }
            return pv;
        }

        private static double ErrorEvaluation(int nn2, Complex[] q2, double ms, double mp, double are, double mre)
        {
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= nn2; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        private static double Cauchy(double[] pt2)
        {
            double xm;
            double f;
            var nn2 = pt2.Length - 1;
            pt2[nn2] = -pt2[nn2];
            var n = nn2;
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / n;
            if (!pt2[n - 1].NearZero())
            {
                xm = -pt2[nn2] / pt2[n - 1];
                if (xm < x) x = xm;
            }
            while (true)
            {
                xm = x * 0.1;
                f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= n - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        private static Complex ComplexDivide(Complex a, Complex b)
        {
            if (b.Real.NearZero() && b.Imaginary.NearZero())
            {
                return new Complex(DblMax, DblMax);
            }
            return a / b;
        }

        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }

    public static class ComplexPolynomialRootFinder5MileStone1
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            var result = new List<Complex>();
            var zeroRoots = CountNumbersZeroRoots(input);
            var copyOfCoeff = new Complex[input.Length - zeroRoots];
            for (var i = 0; i < zeroRoots; i++) result.Add(new Complex(0, 0));
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var degree = copyOfCoeff.Length - 1;
            var _qp = new Complex[degree + 1];
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (copyOfCoeff[0] == new Complex(0, 0)) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var scale = Scale(copyOfCoeff);
            var shift = new Complex(0.70710678, -0.70710678);
            var p = copyOfCoeff.Select(coeff => coeff * scale).ToArray();
            var nn = degree;
            search:
            if (nn <= 1)
            {
                result.Add(ComplexDivide(-p[1], p[0]));
                return result;
            }
            scale = Cauchy(p.Select(ComplexModulus).ToArray());
            for (var cnt1 = 1; cnt1 <= 2; cnt1++)
            {
                var h = NoShift(5, p);
                for (var cnt2 = 1; cnt2 <= 9; cnt2++)
                {
                    shift *= complex90;
                    var (Zero, Converge, qp) = FxShift(10 * cnt2, shift * scale, p, h);
                    if (Converge)
                    {
                        result.Add(Zero.Clone());
                        nn--;
                        p = new Complex[nn + 1];
                        for (var i = 0; i <= nn; i++)
                        {
                            p[i] = qp[i].Clone();
                        }
                        goto search;
                    }
                }
            }
            return result;
        }

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var h = new Complex[nn];
            for (var i = 0; i < nn; i++)
            {
                h[i] = (nn - i) * p[i] / nn;
            }
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(h[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], h[nn - 1]);
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        h[nn - i - 1] = (t * h[nn - i - 2]) + p[nn - i - 1];
                    }
                    h[0] = p[0].Clone();
                }
                else
                {
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        h[nn - i - 1] = h[nn - i - 2];
                    }
                    h[0] = new Complex(0, 0);
                }
            }
            return h;
        }

        private static (Complex Zero, bool Converge, Complex[] qp) FxShift(int l2, Complex s, Complex[] p, Complex[] initH)
        {
            var _qp = new Complex[p.Length];
            var h = initH.Select(elem => elem.Clone()).ToArray();
            var z2 = new Complex(0, 0);
            var pv = PolyEvaluation(s, p, _qp);
            var test = 1;
            var pasd = 0;
            var (hIsZero, t, qh) = Calct(s, pv, h);
            Complex zero;
            bool converge;
            for (var j = 1; j <= l2; j++)
            {
                var ot = t.Clone();
                h = NextH(hIsZero, t, qh, _qp);
                (hIsZero, t, qh) = Calct(s, pv, h);
                z2 = s + t;
                if (!(hIsZero || test != 1 || j == 12))
                {
                    if (ComplexModulus(t - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            (zero, converge, _qp) = Vrshft(10, z2, p, h);
                            if (converge) return (zero, true, _qp);
                            test = 0;
                            pv = PolyEvaluation(s, p, _qp);
                            (hIsZero, t, qh) = Calct(s, pv, h);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            return Vrshft(10, z2, p, h);
        }

        private static (Complex Zero, bool converge, Complex[] qp) Vrshft(int l3, Complex z2, Complex[] p, Complex[] initH)
        {
            var tempH = initH.Select(elem => elem.Clone()).ToArray();
            var nn = p.Length - 1;
            var _qp = new Complex[p.Length];
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            for (int i = 1; i <= l3; i++)
            {
                var pv = PolyEvaluation(s, p, _qp);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(_qp, ms, mp, DblEpsilon, 2.0 * Sqrt(2.0) * DblEpsilon))
                {
                    z2 = s;
                    return (z2, true, _qp);
                }
                (bool HisZero, Complex T, Complex[] qh) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        pv = PolyEvaluation(s, p, _qp);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(s, pv, tempH);
                            tempH = NextH(calcT.HisZero, calcT.T, calcT.qh, _qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }
                    if (mp * 0.1 > omp) return (z2, false, _qp);
                }
                omp = mp;
                _20:
                calcT = Calct(s, pv, tempH);
                tempH = NextH(calcT.HisZero, calcT.T, calcT.qh, _qp);
                calcT = Calct(s, pv, tempH);
                if (!calcT.HisZero)
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false, _qp);
        }

        private static (bool HisZero, Complex T, Complex[] _qh) Calct(Complex s, Complex pv, Complex[] initH)
        {
            var nn = initH.Length - 1;
            var _qh = new Complex[nn + 1];
            var hv = PolyEvaluation(s, initH, _qh);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (bol, t, _qh);
        }

        private static Complex[] NextH(bool bol, Complex t, Complex[] _qh, Complex[] _qp)
        {
            var tempH = new Complex[_qh.Length];
            if (!bol)
            {
                for (var j = 1; j <= _qh.Length - 1; j++)
                {
                    tempH[j] = (t * _qh[j - 1]) + _qp[j];
                }
                tempH[0] = _qp[0];
                return tempH;
            }
            for (var j = 1; j <= _qh.Length - 1; j++)
            {
                tempH[j] = _qh[j - 1];
            }
            tempH[0] = new Complex(0, 0);
            return tempH;
        }

        private static Complex PolyEvaluation(Complex s2, Complex[] p2, Complex[] q2)
        {
            var nn2 = p2.Length - 1;
            q2[0] = p2[0].Clone();
            var pv = q2[0].Clone();
            for (var i = 1; i <= nn2; i++)
            {
                pv = (pv * s2) + p2[i];
                q2[i] = pv.Clone();
            }
            return pv;
        }

        private static double ErrorEvaluation(Complex[] q2, double ms, double mp, double are, double mre)
        {
            var nn2 = q2.Length - 2;
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= nn2; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        private static double Cauchy(double[] pt2)
        {
            double xm;
            double f;
            var nn2 = pt2.Length - 1;
            pt2[nn2] = -pt2[nn2];
            var n = nn2;
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / n;
            if (!pt2[n - 1].NearZero())
            {
                xm = -pt2[nn2] / pt2[n - 1];
                if (xm < x) x = xm;
            }
            while (true)
            {
                xm = x * 0.1;
                f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= n - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        private static Complex ComplexDivide(Complex a, Complex b)
        {
            if (b.Real.NearZero() && b.Imaginary.NearZero())
            {
                return new Complex(DblMax, DblMax);
            }
            return a / b;
        }

        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }

    /// <summary>
    /// Fixed The polyEvaluation function
    /// </summary>
    public static class ComplexPolynomialRootFinder6MileStone2
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            var result = new List<Complex>();
            var zeroRoots = CountNumbersZeroRoots(input);
            var copyOfCoeff = new Complex[input.Length - zeroRoots];
            for (var i = 0; i < zeroRoots; i++) result.Add(new Complex(0, 0));
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var degree = copyOfCoeff.Length - 1;
            var _qp = new Complex[degree + 1];
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (copyOfCoeff[0] == new Complex(0, 0)) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var scale = Scale(copyOfCoeff);
            var shift = new Complex(0.70710678, -0.70710678);
            var p = copyOfCoeff.Select(coeff => coeff * scale).ToArray();
            var nn = degree;
            search:
            if (nn <= 1)
            {
                result.Add(ComplexDivide(-p[1], p[0]));
                return result;
            }
            scale = Cauchy(p.Select(ComplexModulus).ToArray());
            for (var cnt1 = 1; cnt1 <= 2; cnt1++)
            {
                var h = NoShift(5, p);
                for (var cnt2 = 1; cnt2 <= 9; cnt2++)
                {
                    shift *= complex90;
                    var (Zero, Converge, qp) = FxShift(10 * cnt2, shift * scale, p, h);
                    if (Converge)
                    {
                        result.Add(Zero.Clone());
                        nn--;
                        p = new Complex[nn + 1];
                        for (var i = 0; i <= nn; i++)
                        {
                            p[i] = qp[i].Clone();
                        }
                        goto search;
                    }
                }
            }
            return result;
        }

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var h = new Complex[nn];
            for (var i = 0; i < nn; i++)
            {
                h[i] = (nn - i) * p[i] / nn;
            }
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(h[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], h[nn - 1]);
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        h[nn - i - 1] = (t * h[nn - i - 2]) + p[nn - i - 1];
                    }
                    h[0] = p[0].Clone();
                }
                else
                {
                    for (var i = 0; i <= nn - 2; i++)
                    {
                        h[nn - i - 1] = h[nn - i - 2];
                    }
                    h[0] = new Complex(0, 0);
                }
            }
            return h;
        }

        private static (Complex Zero, bool Converge, Complex[] qp) FxShift(int l2, Complex s, Complex[] p, Complex[] initH)
        {
            var h = initH.Select(elem => elem.Clone()).ToArray();
            var z2 = new Complex(0, 0);
            var (pv, qp) = PolyEvaluation(s, p);
            var test = 1;
            var pasd = 0;
            var (hIsZero, t, qh) = Calct(s, pv, h);
            Complex zero;
            bool converge;
            for (var j = 1; j <= l2; j++)
            {
                var ot = t.Clone();
                h = NextH(hIsZero, t, qh, qp);
                (hIsZero, t, qh) = Calct(s, pv, h);
                z2 = s + t;
                if (!(hIsZero || test != 1 || j == 12))
                {
                    if (ComplexModulus(t - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            (zero, converge, qp) = Vrshft(10, z2, p, h);
                            if (converge) return (zero, true, qp);
                            test = 0;
                            (pv, qp) = PolyEvaluation(s, p);
                            (hIsZero, t, qh) = Calct(s, pv, h);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            return Vrshft(10, z2, p, h);
        }

        private static (Complex Zero, bool converge, Complex[] qp) Vrshft(int l3, Complex z2, Complex[] p, Complex[] initH)
        {
            var h = initH.Select(elem => elem.Clone()).ToArray();
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            var qp = new Complex[p.Length];
            Complex pv;
            for (int i = 1; i <= l3; i++)
            {
                (pv, qp) = PolyEvaluation(s, p);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(qp, ms, mp, DblEpsilon, 2.0 * Sqrt(2.0) * DblEpsilon))
                {
                    z2 = s;
                    return (z2, true, qp);
                }
                (bool HisZero, Complex T, Complex[] qh) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        (pv, qp) = PolyEvaluation(s, p);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(s, pv, h);
                            h = NextH(calcT.HisZero, calcT.T, calcT.qh, qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }
                    if (mp * 0.1 > omp) return (z2, false, qp);
                }
                omp = mp;
                _20:
                calcT = Calct(s, pv, h);
                h = NextH(calcT.HisZero, calcT.T, calcT.qh, qp);
                calcT = Calct(s, pv, h);
                if (!calcT.HisZero)
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false, qp);
        }

        private static (bool HisZero, Complex T, Complex[] qh) Calct(Complex s, Complex pv, Complex[] initH)
        {
            var nn = initH.Length - 1;
            var (hv, qh) = PolyEvaluation(s, initH);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (bol, t, qh);
        }

        private static Complex[] NextH(bool bol, Complex t, Complex[] _qh, Complex[] _qp)
        {
            var tempH = new Complex[_qh.Length];
            if (!bol)
            {
                for (var j = 1; j <= _qh.Length - 1; j++)
                {
                    tempH[j] = (t * _qh[j - 1]) + _qp[j];
                }
                tempH[0] = _qp[0];
                return tempH;
            }
            for (var j = 1; j <= _qh.Length - 1; j++)
            {
                tempH[j] = _qh[j - 1];
            }
            tempH[0] = new Complex(0, 0);
            return tempH;
        }

        private static (Complex, Complex[]) PolyEvaluation(Complex s2, Complex[] p2)
        {
            var length = p2.Length;
            var q2 = new Complex[length];
            q2[0] = p2[0].Clone();
            for (var i = 1; i < length; i++)
            {
                q2[i] = (q2[i - 1] * s2) + p2[i];
            }
            return (q2[length - 1], q2);
        }

        private static double ErrorEvaluation(Complex[] q2, double ms, double mp, double are, double mre)
        {
            var nn2 = q2.Length - 2;
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= nn2; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        private static double Cauchy(double[] pt2)
        {
            double xm;
            double f;
            var nn2 = pt2.Length - 1;
            pt2[nn2] = -pt2[nn2];
            var n = nn2;
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / n;
            if (!pt2[n - 1].NearZero())
            {
                xm = -pt2[nn2] / pt2[n - 1];
                if (xm < x) x = xm;
            }
            while (true)
            {
                xm = x * 0.1;
                f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= n - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        private static Complex ComplexDivide(Complex a, Complex b)
        {
            if (b.Real.NearZero() && b.Imaginary.NearZero())
            {
                return new Complex(DblMax, DblMax);
            }
            return a / b;
        }

        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }

        internal static class ComplexPolynomialRootFinderFinal
    {
        private const double DblEpsilon = 2.22044604925031E-16;
        private const double DblMax = 1.79769313486232E+307;
        private const double DblMin = 2.2250738585072E-308;
        private const double DblRadix = 2;

        private static int CountNumbersZeroRoots(Complex[] input)
        {
            var nn = 0;
            while (input[input.Length - 1 - nn] == new Complex(0, 0))
            {
                nn++;
            }
            return nn;
        }

        public static List<Complex> FindRoots(params Complex[] input)
        {
            if (input[0] == 0) throw new Exception("The leading coefficient is zero. No further action taken. Program terminated.");
            var zeroRoots = CountNumbersZeroRoots(input);
            return ((Complex)0).MakeArrayFor(zeroRoots)
                .Concat(
                PartRoots(
                    ScaleCoefficients(input, zeroRoots))).ToList();
        }

        private static Complex[] ScaleCoefficients(Complex[] input, int NumbersOfZeroRoots)
        {
            var copyOfCoeff = new Complex[input.Length - NumbersOfZeroRoots];
            for (var i = 0; i < copyOfCoeff.Length; i++) copyOfCoeff[i] = input[i].Clone();
            var scale = Scale(copyOfCoeff);
            return copyOfCoeff.Select(coeff => coeff * scale).ToArray();
        }

        private static IEnumerable<Complex> PartRoots(Complex[] p)
        {
            var shift = new Complex(0.70710678, -0.70710678);
            var complex90 = new Complex(-0.060756474, -0.99756405);
            if (p.Length <= 2) return new List<Complex> { ComplexDivide(-p[1], p[0]) };
            if (p.Length == 3) return Quadratic(p[0], p[1], p[2]);
            var scale = Cauchy(p.Select(ComplexModulus).ToArray());
            for (var i = 1; i <= 2; i++)
            {
                var h = NoShift(5, p);
                for (var j = 1; j <= 9; j++)
                {
                    shift *= complex90;
                    var (zero, converges, qp) = FxShift(10 * j, shift * scale, p, h);
                    if (converges) return new List<Complex> { zero }.Concat(PartRoots(qp));
                }
            }
            return new List<Complex>();
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

        private static Complex[] NoShift(int l1, Complex[] p)
        {
            var nn = p.Length - 1;
            var h = Derivative(p);
            for (var jj = 1; jj <= l1; jj++)
            {
                if (ComplexModulus(h[nn - 1]) > DblEpsilon * 10 * ComplexModulus(p[nn - 1]))
                {
                    var t = ComplexDivide(-p[nn], h[nn - 1]);
                    for (int i2 = h.Length - 1; i2 >= 1; i2--)
                    {
                        h[i2] = (t * h[i2 - 1]) + p[i2];
                    }
                    h[0] = p[0].Clone();
                }
                else
                {
                    for (int i2 = h.Length - 1; i2 >= 1; i2--)
                    {
                        h[i2] = h[i2 - 1];
                    }
                    h[0] = new Complex(0, 0);
                }
            }
            return h;
        }

        private static (Complex Zero, bool Converge, Complex[] qp) FxShift(int l2, Complex s, Complex[] p, Complex[] initH)
        {
            var h = initH.Select(elem => elem.Clone()).ToArray();
            var z2 = new Complex(0, 0);
            var (pv, qp) = DeflatePolynomial(s, p);
            var test = 1;
            var pasd = 0;
            var (t, qh) = Calct(s, pv, h);
            for (var j = 1; j <= l2; j++)
            {
                var ot = t.Clone();
                h = NextH(t, qh, qp);
                (t, qh) = Calct(s, pv, h);
                z2 = s + t;
                if (!(t == new Complex(0, 0) || test != 1 || j == 12))
                {
                    if (ComplexModulus(t - ot) < 0.5 * ComplexModulus(z2))
                    {
                        if (pasd == 1)
                        {
                            var (zero, converge, qp2) = Vrshft(10, z2, p, h);
                            if (converge) return (zero, true, qp2);
                            test = 0;
                            (pv, qp) = DeflatePolynomial(s, p);
                            (t, qh) = Calct(s, pv, h);
                            continue;
                        }
                        pasd = 1;
                    }
                }
                else
                {
                    pasd = 0;
                }
            }
            return Vrshft(10, z2, p, h);
        }

        private static (Complex Zero, bool converge, Complex[] qp) Vrshft(int l3, Complex z2, Complex[] p, Complex[] initH)
        {
            var h = initH;
            double omp = 0;
            double relstp = 0;
            var b = 0;
            var s = z2.Clone();
            var qp = new Complex[p.Length];
            Complex pv;
            for (int i = 1; i <= l3; i++)
            {
                (pv, qp) = DeflatePolynomial(s, p);
                var mp = ComplexModulus(pv);
                var ms = ComplexModulus(s);
                if (mp <= 20 * ErrorEvaluation(qp, ms, mp))
                {
                    return (s, true, qp);
                }
                (Complex T, Complex[] qh) calcT;
                if (i != 1)
                {
                    if (!(b == 1 || mp < omp || relstp >= 0.05))
                    {
                        var tp = relstp;
                        b = 1;
                        if (relstp < DblEpsilon)
                            tp = DblEpsilon;
                        var r1 = Sqrt(tp);
                        s *= new Complex(1 + r1, r1);
                        (pv, qp) = DeflatePolynomial(s, p);
                        for (var j = 1; j <= 5; j++)
                        {
                            calcT = Calct(s, pv, h);
                            h = NextH(calcT.T, calcT.qh, qp);
                        }
                        omp = DblMax;
                        goto _20;
                    }
                    if (mp * 0.1 > omp) return (z2, false, qp);
                }
                omp = mp;
                _20:
                calcT = Calct(s, pv, h);
                h = NextH(calcT.T, calcT.qh, qp);
                calcT = Calct(s, pv, h);
                if (calcT.T != new Complex(0, 0))
                {
                    relstp = ComplexModulus(calcT.T) / ComplexModulus(s);
                    s += calcT.T;
                }
            }
            return (z2, false, qp);
        }

        private static (Complex t, Complex[] qh) Calct(Complex s, Complex pv, Complex[] initH)
        {
            var nn = initH.Length - 1;
            var (hv, qh) = DeflatePolynomial(s, initH);
            var bol = ComplexModulus(hv) <= DblEpsilon * 10 * ComplexModulus(initH[nn]);
            var t = !bol ? ComplexDivide(-pv, hv) : new Complex(0, 0);
            return (t, qh);
        }

        private static Complex[] NextH(Complex t, Complex[] _qh, Complex[] _qp)
        {
            var h = new Complex[_qp.Length];
            if (t != new Complex(0, 0))
            {
                for (var j = 1; j <= _qh.Length; j++)
                {
                    h[j] = (t * _qh[j - 1]) + _qp[j];
                }
                h[0] = _qp[0];
                return h;
            }
            for (var j = 1; j <= _qh.Length; j++)
            {
                h[j] = _qh[j - 1];
            }
            h[0] = new Complex(0, 0);
            return h;
        }

        private static (Complex, Complex[]) DeflatePolynomial(Complex s2, Complex[] p2)
        {
            var length = p2.Length;
            var q2 = new Complex[length - 1];
            q2[0] = p2[0].Clone();
            for (var i = 1; i < length - 1; i++)
            {
                q2[i] = (q2[i - 1] * s2) + p2[i];
            }
            return ((q2[length - 2] * s2) + p2[length - 1], q2);
        }

        private static double ErrorEvaluation(Complex[] q2, double ms, double mp)
        {
            const double are = DblEpsilon;
            var mre = 2.0 * Sqrt(2.0) * DblEpsilon;
            var e = ComplexModulus(q2[0]) * mre / (are + mre);
            for (var i = 0; i <= q2.Length - 1; i++)
            {
                e = (e * ms) + ComplexModulus(q2[i]);
            }
            return (e * (are + mre)) - (mp * mre);
        }

        private static double Cauchy(double[] pt2)
        {
            var nn2 = pt2.Length - 1;
            pt2[nn2] = -pt2[nn2];
            var x = Exp(Log(-pt2[nn2]) - Log(pt2[0])) / nn2;
            if (!pt2[nn2 - 1].NearZero())
            {
                var xm = -pt2[nn2] / pt2[nn2 - 1];
                if (xm < x) x = xm;
            }
            while (true)
            {
                var xm = x * 0.1;
                var f = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    f = (f * xm) + pt2[i];
                }
                if (f <= 0) break;
                x = xm;
            }
            var dx = x;
            var q = new double[pt2.Length];
            while (Abs(dx / x) > 0.005)
            {
                q[0] = pt2[0];
                for (var i = 1; i <= nn2; i++)
                {
                    q[i] = (q[i - 1] * x) + pt2[i];
                }
                var f = q[nn2];
                var df = q[0];
                for (var i = 1; i <= nn2 - 1; i++)
                {
                    df = (df * x) + q[i];
                }
                dx = f / df;
                x -= dx;
            }
            return x;
        }

        private static double Scale(Complex[] pt)
        {
            var hi = Sqrt(DblMax);
            const double lo = DblMin / DblEpsilon;
            double max = 0;
            var min = DblMax;
            double x;
            for (var i = 0; i <= pt.Length - 1; i++)
            {
                x = ComplexModulus(pt[i]);
                if (x > max) max = x;
                if (!x.NearZero() && x < min) min = x;
            }
            if (min >= lo && max <= hi) return 1;
            x = lo / min;
            double sc;
            if (x <= 1)
            {
                sc = 1 / (Sqrt(max) * Sqrt(min));
            }
            else
            {
                sc = x;
                if (DblMax / sc > max) sc = 1;
            }
            return Pow(DblRadix, ToInt32((Log(sc) / Log(DblRadix)) + 0.5));
        }

        private static Complex ComplexDivide(Complex a, Complex b) =>
            b.Real.NearZero() && b.Imaginary.NearZero() ?
            new Complex(DblMax, DblMax) : a / b;

        private static double ComplexModulus(Complex val)
        {
            var ar = Abs(val.Real);
            var ai = Abs(val.Imaginary);
            if (ar < ai)
            {
                return ai * Sqrt(1.0 + Pow(ar / ai, 2.0));
            }
            if (ar > ai)
            {
                return ar * Sqrt(1.0 + Pow(ai / ar, 2.0));
            }
            return ar * Sqrt(2.0);
        }
    }
}
