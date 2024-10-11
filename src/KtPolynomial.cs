using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using KtExtensions;
using EMDD.KtNumerics;
using static EMDD.KtPolynomials.ComplexPolynomialRootFinder;

namespace EMDD.KtPolynomials
{
    public class Zero : Monomial
    {
        public Zero() : base([0])
        {
        }
    }

    public class One : Monomial
    {
        public One() : base([1])
        {
        }
    }

    public class Monomial : KtPolynomial
    {
        protected internal Monomial(Number[] coeffs) : base(coeffs)
        {
        }

        public static KtPolynomial Create(int degree, Number number)
        {
            if (degree == 0) return Create(number);
            var coeffs = new Number[degree + 1];
            for (int i = 0; i < degree; i++) coeffs[i] = KtComplex.Zero;
            coeffs[degree] = number;
            return new Monomial(coeffs);
        }
    }

    public class KtPolynomial : IEquatable<KtPolynomial>
    {
        private string _variable = "x";
        protected KtPolynomial(params Number[] coeffs)
        {
            Coefficients = coeffs;
        }
        /// <summary>
        /// Create a polynomial
        /// </summary>
        /// <param name="coeffs">array of coefficients starts with the highest exponent</param>
        /// <returns></returns>
        public static KtPolynomial Create(params Number[] coeffs)
        {
            if (coeffs == null || coeffs.Length < 1) return new Zero();
            if (coeffs.Length == 1 && coeffs[0] == 1) return new One();
            var fix = FixCoefficients(coeffs);
            if (fix.Length == 1 && fix[0] == 0) return new Zero();
            if (fix.Length == 1 && fix[0] == 1) return new One();
            return new KtPolynomial(fix);
        }

        /// <summary>
        /// Create a polynomial but replace the coefficient
        /// </summary>
        /// <param name="coeffs">array of coefficients starts with the highest exponent</param>
        /// <returns></returns>
        public static KtPolynomial Create(string variable,params Number[] coeffs)
        {
            if (coeffs == null || coeffs.Length < 1) return new Zero();
            if (coeffs.Length == 1 && coeffs[0] == 1) return new One();
            var fix = FixCoefficients(coeffs);
            if (fix.Length == 1 && fix[0] == 0) return new Zero();
            if (fix.Length == 1 && fix[0] == 1) return new One();
            var poly = new KtPolynomial(fix)
            {
                _variable = variable
            };
            return poly;
        }

        private static Number[] FixCoefficients(Number[] coeffs)
        {
            if (coeffs == null || coeffs.Length < 1) return [KtComplex.Zero];
            var list = coeffs.ToList();
            while (list.Count >= 2 && list[^1] == 0)
            {
                list.RemoveAt(list.Count - 1);
            }
            return list.Select(coeff => coeff.Clone()).ToArray();
        }

        public enum RootMethod
        {
            JenkinsTraub = 1,
            Laguerre = 2,
            Weirstrass = 3
        }

        public Number[] Coefficients { get; }

        public int Degree => Coefficients.Length - 1;

        public static KtPolynomial Integral(KtPolynomial p)
        {
            var buf = new Number[p.Degree + 2];
            buf[0] = KtComplex.Zero;
            for (var i = 1; i < buf.Length; i++) buf[i] = p.Coefficients[i - 1] / i;
            return Create(buf);
        }

        public static double MaxValue(KtPolynomial p, Number[] z) => z.Select(t => p.HornerScheme(t).Magnitude()).Concat(new double[] { 0 }).Max();

        public static KtPolynomial operator -(KtPolynomial p, KtPolynomial q) => p + -q;

        public static KtPolynomial operator -(KtPolynomial p, double num) => p + -num;

        public static KtPolynomial operator -(double num, KtPolynomial p) => num + -p;

        public static KtPolynomial operator -(KtPolynomial p) => Create(p.Coefficients.Select(elem => -elem).ToArray());

        public static bool operator !=(KtPolynomial a, KtPolynomial b) => !(a == b);

        public static KtPolynomial operator *(Number d, KtPolynomial p) => Create(p.Coefficients.Select(elem => d * elem).ToArray());

        public static KtPolynomial operator *(KtPolynomial p, Number d) => Create(p.Coefficients.Select(elem => elem * d).ToArray());

        public static KtPolynomial operator *(KtPolynomial p, KtPolynomial q)
        {
            KtPolynomial r = new Zero();
            for (int i = 0; i <= p.Degree; i++)
            {
                for (int j = 0; j <= q.Degree; j++)
                {
                    r += p.Coefficients[i] * q.Coefficients[j] * Monomial.Create(i + j, 1);
                }
            }
            return r;
        }

        public static KtPolynomial operator /(KtPolynomial p, Number d) => p * d.Inverse();

        public static (KtPolynomial Numerator, KtPolynomial Denominator) operator /(KtPolynomial p1, KtPolynomial p2) => (p1, p2) switch
        {
            (null, _) => (new Zero(), new One()),
            (Zero zero, _) => (zero, new One()),
            (_, null) => (new One(), new Zero()),
            (_, Zero zero) => (new One(), zero),
            (One u, _) => (u, p2.Clone()),
            (_, One u) => (p1.Clone(), u),
            (_, { Degree: 0 }) => (p1 / p2.Coefficients[0], new One()),
            _ => Try(p1, p2)
        };

        private static (KtPolynomial Numerator, KtPolynomial Denominator) Try(KtPolynomial p1, KtPolynomial p2)
        {
            var commonRoots = p1.Roots().Intersect(p2.Roots()).ToArray();
            var numerator = commonRoots.Aggregate(p1, (current, number) => current.Deflate(-number));
            var denominator = commonRoots.Aggregate(p2, (poly, root) => poly.Deflate(-root));
            var rationalizer = NumberHelper.Gcd(numerator.LeadingCoefficient(), denominator.LeadingCoefficient());
            return (numerator / rationalizer, denominator / rationalizer);
        }

        public static KtPolynomial operator ^(KtPolynomial p, uint k)
        {
            if (k == 0) return new One();
            if (k == 1) return p;
            return p * (p ^ (k - 1));
        }

        public static KtPolynomial operator +(KtPolynomial p, KtPolynomial q)
        {
            var degree = Math.Max(p.Degree, q.Degree);
            var coeffs = new Number[degree + 1];
            for (var i = 0; i <= degree; i++)
            {
                if (i > p.Degree) coeffs[i] = q.Coefficients[i];
                else if (i > q.Degree) coeffs[i] = p.Coefficients[i];
                else coeffs[i] = p.Coefficients[i] + q.Coefficients[i];
            }
            return Create(coeffs);
        }

        public static KtPolynomial operator +(KtPolynomial p, double num) => p + Monomial.Create(0, num);

        public static KtPolynomial operator +(double num, KtPolynomial p) => Monomial.Create(0, num) + p;

        public static bool operator ==(KtPolynomial a, KtPolynomial b)
        {
            if (ReferenceEquals(a, b)) return true;
            if (a is null || b is null) return false;
            return a.Degree == b.Degree && a.Coefficients.SequenceEqual(b.Coefficients);
        }

        public KtPolynomial Clone() => Create(Coefficients.Select(elem => elem).ToArray());

        public KtPolynomial Deflate(Number val)
        {
            var copy = Coefficients.Reverse().ToArray();
            var result = new Number[Degree];
            result[0] = copy[0];
            for (var i = 1; i < Degree; i++)
                result[i] = (result[i - 1] * -val) + copy[i];
            return Create(result.Reverse().ToArray());
        }

        public KtPolynomial Derivative()
        {
            if (Degree <= 0) return new Zero();
            var buf = new Number[Degree];
            for (var i = 0; i < buf.Length; i++) buf[i] = (i + 1) * Coefficients[i + 1];
            return Create(buf);
        }

        public KtPolynomial Derivative(int times)
        {
            var clone = Clone();
            for (var i = 0; i < times; i++)
                clone = clone.Derivative();
            return clone;
        }

        public Number Differentiate(Number x) => Derivative().HornerScheme(x);

        public Number Differentiate(Number x, int degree) => Derivative(degree).HornerScheme(x);

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(this, obj)) return true;
            if (this is null || obj is null) return false;
            return obj is KtPolynomial poly && Equals(poly);
        }

        public bool Equals(KtPolynomial other)
        {
            if (ReferenceEquals(this, other)) return true;
            if (this is null || other is null) return false;
            return this == other;
        }

        public override int GetHashCode() => ToString().GetHashCode();

        public Number HornerScheme(Number x)
        {
            var buf = Coefficients[Degree];
            for (var i = Degree - 1; i >= 0; i--)
                buf = Coefficients[i] + (x * buf);
            return buf;
        }

        public Number Integrate(Number a, Number b)
        {
            var buf = new Number[Degree + 2];
            buf[0] = KtComplex.Zero;
            for (var i = 1; i < buf.Length; i++)
                buf[i] = Coefficients[i - 1] / i;
            var p = Create(buf);
            return p.HornerScheme(b) - p.HornerScheme(a);
        }

        public bool IsAZero(Number num, int tolerance = 15)
        {
            if (num == null) return false;
            if (num.Magnitude().InvalidNumber()) return false;
            return HornerScheme(num).NearZero(tolerance);
        }

        public Number LeadingCoefficient()
        {
            if (Degree >= 0) return Coefficients[Degree];
            throw new OverflowException("Invalid numbers of coefficients");
        }

        public KtPolynomial Normalize() => Create(DivideEachCoeffByFirstTerm());

        private Number[] DivideEachCoeffByFirstTerm() =>
            Coefficients.Select(elem => elem / Coefficients[Degree]).ToArray();

        public Number[] Roots(RootMethod method = RootMethod.JenkinsTraub)
        {
            if (Degree < 1) return [];
            var roots = GetInitialRootByMethod(method);
            FixRootRoundingErrors(roots);
            return roots;
        }

        private void FixRootRoundingErrors(Number[] roots)
        {
            for (var i = 0; i < roots.Length; i++)
            {
                roots[i] = roots[i].Round(15);
                roots[i] = RoundRootsBaseOnXzPlaneApproach(roots[i]);
                roots[i] = roots[i].SmartRound();
            }
        }

        private Number[] GetInitialRootByMethod(RootMethod method) => method switch
        {
            RootMethod.JenkinsTraub => RootsJenkinsTraub(),
            RootMethod.Laguerre => RootsLaguerre(),
            RootMethod.Weirstrass => RootsWeirstrass(),
            _ => throw new ArgumentOutOfRangeException(nameof(method), method, "Root Method is not implemented yet"),
        };

        public Number[] RootsJenkinsTraub() => FindRoots(Coefficients.Select(coeff => (Complex)coeff).Reverse().ToArray()).Select(root => (Number)root).ToArray();

        public Number[] RootsLaguerre()
        {
            const int tolerance = 8;
            var roots = new List<Number>();
            var curPol = Clone();
            while (roots.Count < Degree)
            {
                Number initial = 0;
                while (true)
                {
                    if (IsAZero(initial, tolerance))
                    {
                        initial = initial.Round(tolerance);
                        curPol = curPol.Deflate(-initial);
                        roots.Add(initial);
                        break;
                    }
                    var degree = curPol.Degree;
                    var evaluation = curPol.HornerScheme(initial);
                    var g = curPol.Differentiate(initial) / evaluation;
                    var h = (g * g) - (2 * curPol.Differentiate(initial, 2) / evaluation);
                    var determ = ((degree - 1) * ((degree * h) - (g * g))).Sqrt();
                    var denom1 = g - determ;
                    var denom2 = g + determ;
                    var a = degree / Number.Max(denom1, denom2);
                    initial -= a;
                }
            }
            return [.. roots];
        }

        public Number[] RootsWeirstrass()
        {
            const int toleranceLevel = 15;
            var tolerance = Math.Pow(10, -toleranceLevel);
            const int maxIterations = 100;
            var normalized = Normalize();
            var degree = normalized.Degree;
            var z = new Number[degree];
            var w = new Number[degree];
            for (var k = 0; k < normalized.Degree; k++)
                z[k] = KtComplex.CreateComplexNumberFromPolar(1, 360.0 * k / normalized.Degree);
            for (var iter = 0; iter < maxIterations && MaxValue(normalized, z) > tolerance; iter++)
            {
                for (var i = 0; i < 10; i++)
                {
                    for (var k = 0; k < degree; k++)
                    {
                        var evaluation = normalized.HornerScheme(z[k]);
                        var weierNull = WeierNull(z, k);
                        w[k] = evaluation / weierNull;
                    }
                    for (var k = 0; k < degree; k++) z[k] -= w[k];
                }
            }

            for (int k = 0; k < normalized.Degree; k++)
            {
                z[k] = z[k].Round();
            }

            return z;
        }

        public override string ToString() => TurnCoeffToString(Coefficients);

        internal string TurnCoeffToString(Number[] coeffsSource)
        {
            if (coeffsSource == null) return "";
            var coeffs = coeffsSource.Select(elem => elem.Clone()).Reverse().ToArray();
            var degree = coeffs.Length - 1;
            var builder = new StringBuilder();
            for (var i = 0; i <= degree; i++)
            {
                var coefficient = coeffs[i];
                if (coefficient == 0) continue;
                builder.Append(ConvertCoefficientToString(coefficient, i, degree));
                if (i < degree) builder.Append(i == degree - 1 ? _variable : $"{_variable}^{degree - i}");
            }
            return builder.ToString();
        }

        internal static Number WeierNull(Number[] z, int k)
        {
            if (k < 0 || k >= z.Length) throw new ArgumentOutOfRangeException(nameof(k));
            Number buf = 1;
            for (var j = 0; j < z.Length; j++)
                if (j != k) buf *= z[k] - z[j];
            return buf;
        }

        private static string ConvertCoefficientToString(Number coefficient, int index, int degree)
        {
            if (index == 0 && degree == 0) return coefficient;
            if (coefficient is KtComplex) return $"+({coefficient})";
            var sign = coefficient < 0 ? "-" : index > 0 ? "+" : "";
            var abscoeff = coefficient.Abs();
            if (index == degree) return $"{sign}{abscoeff}";
            if (abscoeff == 1) return sign;
            return $"{sign}{abscoeff}";
        }

        private Number RoundRootsBaseOnXzPlaneApproach(Number root)
        {
            var evaluation = new double[] { 10000000, 10000000 };
            for (var j = 1; j < 15; j++)
            {
                evaluation[0] = HornerScheme(root.Round(15 - j)).Magnitude();
                if (evaluation[0].NearZero(7))
                {
                    return root.Round(15 - j);
                }
                if (evaluation[0] > evaluation[1])
                {
                    return root.Round(16 - j);
                }
                if (j == 15)
                {
                    return root.Round(15 - j);
                }
                evaluation[1] = evaluation[0];
            }
            return root;
        }
    }
}
