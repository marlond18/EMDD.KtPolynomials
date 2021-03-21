using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using KtExtensions;
using EMDD.KtNumerics;
using EMDD.KtPolynomials;
using static System.Math;
using static EMDD.KtPolynomials.ComplexPolynomialRootFinder;

namespace KtPolynomialTest
{
    public static class Generator
    {
        private static readonly Random _rand = new(DateTime.Now.Millisecond);

        public static int RandomNumber(int lower, int upper)
        {
            return _rand.Next(lower, upper);
        }

        public static Complex Round(Complex number, int precision)
        {
            return new Complex(Math.Round(number.Real, precision), Math.Round(number.Imaginary, precision));
        }

        public static void DebugArray<T>(List<T> val)
        {
            var stringToReturn = val.Aggregate("",
                (current, ktNumericBase) => $"{current}{(current?.Length == 0 ? "" : ",")}{ktNumericBase}");
            Debug.Print(stringToReturn);
        }

        public static void DebugArray<T>(T[] val)
        {
            var stringToReturn = val.Aggregate("",
                (current, ktNumericBase) => $"{current}{(current?.Length == 0 ? "" : ",")}{ktNumericBase}");
            Debug.Print(stringToReturn);
        }
    }

    [TestClass]
    public class KtPioynomialTest2
    {
        [TestMethod]
        public void TestJenkinsTraub2()
        {
            var coeff = new Complex[] { 1, -3, -2, 6 };
            var expectedRoots = new double[] { Round(Sqrt(2), 9), -Round(Sqrt(2), 9), 3 };
            var root = FindRoots(coeff);
            var realroots = root.ConvertAll(elem => Generator.Round(elem ,9).Real );
            Assert.AreEqual(expectedRoots.Length, realroots.Count);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        private static bool SequenceEqual<T>(T[] expected, T[] actual)
        {
            var e = expected.ToList();
            e.Sort();
            var a = actual.ToList();
            a.Sort();
            return a.SequenceEqual(e) && e.SequenceEqual(a);
        }

        private static bool ContentEqual<T>(T[] expected, T[] actual)
        {
            return expected.All(actual.Contains) && actual.All(expected.Contains);
        }

        [TestMethod]
        public void TestJenkinsTraub5()
        {
            var coeff = new Complex[] { 5, -20, 5, 50, -20, -40 };
            var expectedRoots = new double[] { -1, -1, 2, 2, 2.0 };
            var root = FindRoots(coeff);
            var realroots = root.Select(elem => Generator.Round(elem,7).Real).ToArray();
            Assert.AreEqual(expectedRoots.Length, realroots.Length);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub6()
        {
            var coeff = new Complex[] { 1, -4, -18, 108, -135, 0, 0, 0, 0 };
            var expectedRoots = new double[] { 0, 0, 0, 0, -5.0, 3, 3, 3 };
            var root = FindRoots(coeff);
            var realroots = root.Select(elem => Generator.Round(elem,3).Real).ToArray();
            Assert.AreEqual(expectedRoots.Length, realroots.Length);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub7()
        {
            var coeff = new Complex[] { 1, 10, 27, 0, -57, -30, 29, 20 };
            var expectedRoots = new double[] { -1, -1, -1, 1, 1, -5, -4.0 };
            var root = FindRoots(coeff);
            var realroots = root.Select(elem => Generator.Round(elem,7).Real).ToArray();
            Assert.AreEqual(expectedRoots.Length, realroots.Length);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub8()
        {
            var coeff = new Complex[] { 1, -30, -9744, 296918, -558657, -28305288 };
            var expectedRoots = new double[] { -8, 19, 19, 99, -99.0 };
            var root = FindRoots(coeff);
            var realroots = root.Select(elem => Generator.Round(elem,7).Real).ToArray();
            Assert.AreEqual(expectedRoots.Length, realroots.Length);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub()
        {
            var coeff = new Complex[] { 1, -1.5, -5.5, 3 };
            var expectedRoots = new double[] { 0.5, -2, 3 };
            var root = FindRoots(coeff);
            var realroots = root.Select(elem =>Generator.Round(elem, 9).Real).ToArray();
            Assert.AreEqual(expectedRoots.Length, realroots.Length);
            Assert.IsTrue(SequenceEqual(realroots, expectedRoots));
            Assert.IsTrue(ContentEqual(realroots, expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraubWithComplexRoots()
        {
            var coeff = new Complex[] { 1, -2, 2 };
            var expectedRoots = new Complex[] { new Complex(1, -1), new Complex(1, 1) };
            var root = FindRoots(coeff);
            var realroots = root.ConvertAll(elem => Generator.Round(elem,9));
            Assert.AreEqual(expectedRoots.Length, realroots.Count);
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraubRandomNumber()
        {
            for (var i = 0; i < 200; i++)
            {
                var count = Generator.RandomNumber(2, 10);
                var coeffs = new double[count];
                coeffs = coeffs.Select(_ => (double)Generator.RandomNumber(-100, 100)).ToArray();
                var poly = coeffs.Aggregate((KtPolynomial)new One(),
                    (current, ktNumericBase) => current * KtPolynomial.Create(-ktNumericBase,1));
                var roots = poly.Roots();
                var fixroots = roots.Select(elem => elem.Round(7).ToComplex().Real).ToArray();
                Assert.AreEqual(coeffs.Length, fixroots.Length);
                if (!SequenceEqual(coeffs, fixroots))
                {
                    Debug.Print("Polynomial");
                    Debug.Print(poly.ToString());
                    Debug.Print("expected");
                    Generator.DebugArray(coeffs);
                    Debug.Print("actual");
                    Generator.DebugArray(fixroots);
                    Assert.Fail();
                }
                Assert.IsTrue(SequenceEqual(coeffs, fixroots));
                Assert.IsTrue(ContentEqual(coeffs, fixroots));
            }
        }
    }

    [TestClass]
    public class KtPolyTraubJenkins
    {
        [TestMethod]
        public void TestJenkinsTraub2()
        {
            var coeff = new Complex[] { 1, -3, -2, 6 };
            Number[] expectedRoots = { Round(Sqrt(2), 9), -Round(Sqrt(2), 9), 3 };
            var roots = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).Roots();
            Assert.AreEqual(expectedRoots.Length, roots.Length);
            Assert.IsTrue(roots.Contains(expectedRoots[0]));
            Assert.IsTrue(roots.All(expectedRoots.Contains));
            Assert.IsTrue(expectedRoots.All(roots.Contains));
            Assert.IsTrue(ContentEqual(roots.ToArray(), expectedRoots));
        }

        private static bool SequenceEqual<T>(T[] expected, T[] actual)
        {
            var e = expected.ToList();
            e.Sort();
            var a = actual.ToList();
            a.Sort();
            return a.SequenceEqual(e) && e.SequenceEqual(a);
        }

        private static bool ContentEqual<T>(T[] expected, T[] actual)
        {
            return expected.All(actual.Contains) && actual.All(expected.Contains);
        }

        [TestMethod]
        public void TestJenkinsTraub5()
        {
            var coeff = new Complex[] { 5, -20, 5, 50, -20, -40 };
            var expectedRoots = new Number[] { -1, -1, 2, 2, 2.0 };
            var roots = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).Roots();
            Assert.AreEqual(expectedRoots.Length, roots.Length);
            Assert.IsTrue(SequenceEqual(roots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(roots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestRoot5746()
        {
            var coeff = new Number[]{ 6, 4, -7, 5 };
            var expectedRoots = new Number[] {(1,1),(1,-1),-3.0/5 };
            var roots = KtPolynomial.Create(coeff).Roots();
            Assert.AreEqual(expectedRoots.Length, roots.Length);
            Assert.IsTrue(ContentEqual(roots.ToArray(), expectedRoots));
            Assert.IsTrue(SequenceEqual(roots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void SmartRounding()
        {
            var d = new KtComplex(2.00000003, -1.34000000005E-12);
            Assert.AreEqual(d.SmartRound(), (2, -1.34E-12));
        }

        [TestMethod]
        public void TestForIenumTake()
        {
            var array = new[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            var array2 = new[] { 1, 2, 3, 4 };
            Assert.IsTrue(array.Take(4).SequenceEqual(array2));
        }

        private class Duck
        {
            public readonly int _y;

            public Duck(int _y)
            {
                this._y = _y;
            }
        }

        [TestMethod]
        public void TestIEAggregate()
        {
            var array = new[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            var d = array.Aggregate(new Duck(0), (_, f) =>new Duck(f));
        }

        [TestMethod]
        private void TestTupleEquality<T,TD>((T,TD) d, (T,TD) e)
        {
            Assert.AreEqual(d,e);
            Assert.IsTrue(Equals(d, e) );
            Assert.IsTrue(d.Equals(e));
            Assert.AreEqual(d.GetHashCode(),e.GetHashCode());
        }

        [TestMethod]
        public void MultipleTupleEqualityTest()
        {
            TestTupleEquality(((Number)2, (Number)4), (2, 4));
            TestTupleEquality((1.2, 2.4), (1.2, 2.4));
        }

        [TestMethod]
        public void TestJenkinsTraub6()
        {
            var coeff = new Complex[] { 1, -4, -18, 108, -135, 0, 0, 0, 0 };
            var expectedRoots = new Number[] { 0, 0, 0, 0, -5.0, 3, 3, 3 };
            var root = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).RootsJenkinsTraub();
            var realroots = root.Select(elem => elem.Round(3)).ToList();
            Assert.AreEqual(expectedRoots.Length, realroots.Count);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub7()
        {
            var coeff = new Complex[] { 1, 10, 27, 0, -57, -30, 29, 20 };
            var expectedRoots = new Number[] { -1, -1, -1, 1, 1, -5, -4.0 };
            var roots = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).Roots();
            Assert.AreEqual(expectedRoots.Length, roots.Length);
            Assert.IsTrue(SequenceEqual(roots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(roots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub8()
        {
            var coeff = new Complex[] { 1, -30, -9744, 296918, -558657, -28305288 };
            var expectedRoots = new Number[] { -8, 19, 19, 99, -99.0 };
            var root = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).RootsJenkinsTraub();
            var realroots = root.Select(elem => elem.Round(7)).ToList();
            Assert.AreEqual(expectedRoots.Length, realroots.Count);
            Assert.IsTrue(SequenceEqual(realroots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(realroots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraub()
        {
            var coeff = new Complex[] { 1, -1.5, -5.5, 3 };
            var expectedRoots = new Number[] { 0.5, -2, 3 };
            var roots = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).Roots();
            Assert.AreEqual(expectedRoots.Length, roots.Length);
            Assert.IsTrue(SequenceEqual(roots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(roots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraubWithComplexRoots()
        {
            var coeff = new Complex[] { 1, -2, 2 };
            var expectedRoots = new[] { (Number)new Complex(1, -1), (Number)new Complex(1, 1) };
            var roots = KtPolynomial.Create(coeff.Select(elem => (Number)elem).Reverse().ToArray()).Roots();
            Assert.AreEqual(expectedRoots.Length, roots.Length);
            Assert.IsTrue(SequenceEqual(roots.ToArray(), expectedRoots));
            Assert.IsTrue(ContentEqual(roots.ToArray(), expectedRoots));
        }

        [TestMethod]
        public void TestJenkinsTraubRandomNumber()
        {
            var fails = 0;
            for (var i = 0; i < 200; i++)
            {
                var count = Generator.RandomNumber(2, 20);
                var coeffs = new Number[count];
                coeffs = coeffs.Select(_ => (Number)Generator.RandomNumber(-40, 100)).ToArray();
                var poly = coeffs.Aggregate((KtPolynomial)new One(),
                    (current, ktNumericBase) => current * KtPolynomial.Create(-ktNumericBase,1));
                var roots = poly.Roots();
                Assert.AreEqual(coeffs.Length, roots.Length);
                if (!SequenceEqual(coeffs, roots))
                {
                    fails++;
                    Debug.Print("Polynomial");
                    Generator.DebugArray(poly.Coefficients.Reverse().ToArray());
                    Debug.Print("expected");
                    Generator.DebugArray(coeffs);
                    Debug.Print("actual");
                    Generator.DebugArray(roots);
                }
                else
                {
                    Debug.Print("Passed");
                    Generator.DebugArray(roots);
                }
                //Assert.IsTrue(SequenceEqual(coeffs, roots));
                //Assert.IsTrue(ContentEqual(coeffs, roots));
            }
            Debug.Print($"{fails}");
        }
    }

    [TestClass]
    public class DividingPolynomials
    {
        [TestMethod]
        public void Division1()
        {
            var poly1 = KtPolynomial.Create(4,4,1);
            var poly2 = KtPolynomial.Create(2,1);
            var (Numerator, Denominator) = poly1 / poly2;
            Assert.AreEqual(Numerator,poly2);
            Assert.AreEqual(Denominator, new One());
        }

        [TestMethod]
        public void Division2()
        {
            var poly1 = KtPolynomial.Create(4, -3,-4,3);
            var poly2 = KtPolynomial.Create(-48, 88, 53, -77, 6);
            var (Numerator, Denominator) = poly1 / poly2;
            Assert.AreEqual(Numerator, KtPolynomial.Create(-1,1));
            Assert.AreEqual(Denominator, KtPolynomial.Create(12,-25,2));
        }

        [TestMethod]
        public void Division3()
        {
            var poly1 = KtPolynomial.Create(-10, 1,0, 4, 1);
            var poly2 = KtPolynomial.Create(-5, 3, 1);
            var (Numerator, Denominator) = poly1 / poly2;
            Assert.AreEqual(Numerator,   KtPolynomial.Create(2,1,1));
            Assert.AreEqual(Denominator, KtPolynomial.Create(1));
        }

        [TestMethod]
        public void DivisionOfZero()
        {
            var poly1 = KtPolynomial.Create(0);
            var poly2 = KtPolynomial.Create(-5, 3, 1);
            var (Numerator, Denominator) = poly1 / poly2;
            Assert.AreEqual(Numerator,  KtPolynomial.Create(0));
            Assert.AreEqual(Denominator,KtPolynomial.Create(1));
        }

        [TestMethod]
        public void Division4()
        {
            var poly2 = KtPolynomial.Create(-1,0,-7,0,10);
            var poly1 = KtPolynomial.Create(3,-1,1);
            var (Numerator, Denominator) = poly1 / poly2;
            Assert.AreEqual(Denominator, poly2);
            Assert.AreEqual(Numerator, poly1);
        }

        [TestMethod]
        public void Division5()
            {
            var poly2 = KtPolynomial.Create(0,0,0,16,0,0,36,0,4,18,0,6);
            var poly1 = KtPolynomial.Create(4,0,0,3,0,1);
            var (Numerator, Denominator) = poly1 / poly2;
            Assert.AreEqual(Numerator,   KtPolynomial.Create(1));
            Assert.AreEqual(Denominator, KtPolynomial.Create(0,0,0,4,0,0,6));
        }
    }
}
