using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.WindowsForms;

namespace Lab2
{
    public delegate double someFunction(double x);

    public class CubicSpline
    {
        public class Polinom
        {
            public double X;
            public double a, b, c, d;
            public Polinom(double x, double a, double b, double c, double d)
            {
                this.X = x;
                this.a = a;
                this.b = b;
                this.c = c;
                this.d = d;

                //Console.WriteLine($"\n a = {a}, b = {b}, c = {c}, d = {d}");

            }

            public double Calculate(double x)
            {
                return a + b * (x - this.X) + c * Math.Pow((x - this.X), 2) / 2 + d * Math.Pow((x - this.X), 3) / 6;
            }
        }

        public Polinom[] Polinoms;

        public double[] ControlPoints;

        public double Calculate(double x)
        {
            //Console.WriteLine($"Calculating for {x}");

            int polinomIndex = -1;

            for (int i = 0; i < ControlPoints.Length - 1; i++)
            {
                if (x >= ControlPoints[i] && x <= ControlPoints[i + 1])
                {
                    polinomIndex = i;

                    break;
                }
            }

            if (polinomIndex == -1)
            {
                //throw new ArgumentException($"Argument {x} is ot of spline range");

                polinomIndex = ControlPoints.Length - 2;
            }

            //Console.WriteLine($"Used peace {polinomIndex + 1}");

            double result = Polinoms[polinomIndex].Calculate(x);

            //Console.WriteLine($"Got result {result}");

            return result;

        }


        public CubicSpline(double[] controlPoints, double[] controlResults, double b0, double bn)
        {
            double[] a = new double[controlPoints.Length];
            double[] b = new double[controlPoints.Length];
            double[] c = new double[controlPoints.Length];
            double[] d = new double[controlPoints.Length];

            double[] h = new double[controlPoints.Length];

            for (int i = 0; i < controlPoints.Length; i++)
            {
                //Console.WriteLine($"x = {controlPoints[i]}, y = {controlResults[i]}");
            }

            ControlPoints = controlPoints;

            for (int i = 0; i < a.Length; i++)
            {
                a[i] = controlResults[i];
            }

            for (int i = 1; i < h.Length; i++)
            {
                h[i] = ControlPoints[i] - controlPoints[i - 1];
            }

            double[] someC = GetCValues(h, controlResults);

            for (int i = 1; i < c.Length - 1; i++)
            {
                c[i] = someC[i - 1];
            }

            c[0] = (6 / h[1]) * (b0 - ((a[1] - a[0]) / h[1])) - 2 * c[1];

            c[^1] = ((6 / h[^1]) * (bn - ((a[^1] - a[^2]) / h[^1])) - c[^2]) / 2;

            for (int i = 1; i < d.Length; i++)
            {
                d[i] = (c[i] - c[i - 1]) / h[i];
            }

            b[^1] = bn;

            for (int i = 1; i < b.Length - 1; i++)
            {
                b[i] = ((a[i] - a[i - 1]) / h[i]) + (2 * c[i] + c[i - 1]) * h[i] / 6;
            }

            Polinoms = new Polinom[controlPoints.Length - 1];

            for (int i = 0; i < a.Length - 1; i++)
            {
                Polinoms[i] = new Polinom(controlPoints[i + 1], a[i + 1], b[i + 1], c[i + 1], d[i + 1]);
            }

        }

        /// <summary>
        /// Gets gamma coefs ('c' coefs in code code) by solving tridiagonal matrix.
        /// </summary>
        /// <param name="h"></param>
        /// <param name="controlResults"></param>
        /// <returns></returns>
        private double[] GetCValues(double[] h, double[] controlResults)
        {
            double[] F = new double[controlResults.Length - 2];

            for (int i = 0; i < F.Length; i++)
            {
                F[i] = 6 * ((((controlResults[i + 1 + 1] - controlResults[i + 1]) / h[i + 1 + 1]) - (controlResults[i + 1] - controlResults[i - 1 + 1]) / h[i + 1])) / (h[i + 1] + h[i + 1 + 1]);
            }


            double[,] matrix = new double[F.Length, F.Length];

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, i] = 2;
            }

            for (int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                matrix[i, i + 1] = h[i + 2] / (h[i + 1] + h[i + 1 + 1]);
            }

            for (int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                matrix[i + 1, i] = h[i + 2] / (h[i + 2] + h[i + 3]);
            }

            return SolveTriDiagnlMatrix(matrix, F);
        }

        private double[] SolveTriDiagnlMatrix(double[,] matrix, double[] f)
        {
            matrix = (double[,])matrix.Clone();
            f = (double[])f.Clone();

            for (int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                double a = (matrix[i + 1, i] / matrix[i, i]);

                matrix[i + 1, i + 1] += -a * matrix[i, i + 1];

                f[i + 1] += -a * f[i];
            }

            double[] result = new double[f.Length];


            result[^1] = f[^1] / matrix[matrix.GetLength(0) - 1, matrix.GetLength(0) - 1];

            int mSize = matrix.GetLength(0);

            for (int i = f.Length - 2; i >= 0; i--)
            {
                //result[^(i + 1)] = (f[i] - matrix[matrix.GetLength(0) - 1 - i, matrix.GetLength(0) - 1 - i + 1] * result[^i]) / matrix[matrix.GetLength(0) - 1 - i, matrix.GetLength(0) - 1 - i ];

                result[i] = (f[i] - matrix[i, i + 1] * result[i + 1]) / (matrix[i, i]);
            }

            return result;

        }
    }

    class Program
    {
        static void Main(string[] args)
        {
            DoAnalytic(func1, func1Devir, 0, 3,
                8);

            //DoConvergenceAnalytic(func1, func1Devir, 0, 3);
        }

        /// <summary>
        /// Builds new spline and draws its plot.
        /// </summary>
        /// <param name="funcToInterpolate"></param>
        /// <param name="funcDev"></param>
        /// <param name="xmin"></param>
        /// <param name="xmax"></param>
        /// <param name="pointsNumber"></param>
        public static void DoAnalytic(Func<double, double> funcToInterpolate, Func<double, double> funcDev, double xmin, double xmax, int pointsNumber)
        {
            double[] xPoints = GetControlPoints(xmin, xmax, pointsNumber).ToArray();
            double[] yPoints = xPoints.Select((d, i) => funcToInterpolate(d)).ToArray();

            //double[] xPoints = new[] {1.0, 2, 3, 4, 5, 6};
            //double[] yPoints = new[] { 1.0002, 1.0341, 0.6, 0.40105, 0.1, 0.23975 };

            CubicSpline spline = new CubicSpline(xPoints, yPoints, funcDev(xmin), funcDev(xmax));

            double error = GetError(spline.Calculate, funcToInterpolate, xmin, xmax);

            Console.WriteLine($"Error = {error}");

            ShowGraph(spline.Calculate, xmin, xmax, 0.01, "Spline", ToControlPoints(xPoints, yPoints), funcToInterpolate);

            ShowGraph(funcOst(funcToInterpolate, spline.Calculate), xmin, xmax, 0.001, "f(x) - s(x)");
        }

        /// <summary>
        /// Draws spline coverage histogram
        /// </summary>
        /// <param name="funcToInterpolate"></param>
        /// <param name="funcDev"></param>
        /// <param name="xmin"></param>
        /// <param name="xmax"></param>
        public static void DoConvergenceAnalytic(Func<double, double> funcToInterpolate, Func<double, double> funcDev, double xmin, double xmax)
        {
            List<DataPoint> points = new List<DataPoint>();

            for (int pointsCount = 3; pointsCount < 1000; pointsCount+=10)
            {

                double[] xPoints = GetControlPoints(xmin, xmax, pointsCount).ToArray();
                double[] yPoints = xPoints.Select((d, i) => funcToInterpolate(d)).ToArray();

                CubicSpline spline = new CubicSpline(xPoints, yPoints, funcDev(xmin), funcDev(xmax));

                double error = GetError(spline.Calculate, funcToInterpolate, xmin, xmax);

               points.Add(new DataPoint(pointsCount,( error)));
            }

            ShowGraph(points);
        }

        public static List<double> GetControlPoints(double xMin, double xMax, int number)
        {
            var result = new List<double>();

            double step = (xMax - xMin) / (number - 1);

            for (int i = 0; i < number; i++)
            {
                result.Add(xMin + i * step);
            }

            return result;
        }

        public static double GetError(Func<double,double> realFunc, Func<double, double> interpolatedFunc, double xmin, double xmax)
        {
            List<double> points = GetControlPoints(xmin, xmax, 1013);

            Func<double, double> diffFunc = funcOst(realFunc, interpolatedFunc);

            return points.Select(point =>Math.Abs(diffFunc(point))).Max();
        }

        public static List<double> GetControlPoints(double xMin, double xMax, double step)
        {
            var result = new List<double>();

            for (double i = 0; i <= xMax - xMin; i += step)
            {
                result.Add(xMin + i);
            }

            return result;
        }


        public static double func1(double x)
        {
            return (x * x) * Math.Sin(2 * x) / (2 + Math.Cos(5 * x));
        }

        public static double func2(double x)
        {
            return 1 / (1 + 25 * (x - 2) * (x - 2));
        }

        public static double func1Devir(double x)
        {
            return (((x * (2 * x * Math.Cos(2 * x) * (2 + Math.Cos(5 * x)) + Math.Sin(2 * x) * (5 * x * Math.Sin(5 * x) + 2 * Math.Cos(5 * x) + 4))) / Math.Pow(2 + Math.Cos(5 * x), 2)));
        }

        public static double func2Devir(double x)
        {
            return (50 * x - 99) / Math.Pow((25 * (x - 2) * (x - 2) + x), 2);
        }

        public static Func<double, double> funcOst(Func<double, double> f1, Func<double, double> f2)
        {
            return new Func<double, double>(a => f1(a) - f2(a));
        }

        public static double FuncTest(double x)
        {
            return Math.Sin(4 * x);
        }


        public static List<DataPoint> ToControlPoints(double[] x, double[] y)
        {
            var points = new List<DataPoint>();

            for (int i = 0; i < x.Length; i++)
            {
                points.Add(new DataPoint(x[i], y[i]));
            }


            return points;
        }

        public static void ShowGraph(List<DataPoint> dataPoints)
        {

            Form graphForm = new Form();

            PlotView plot = new PlotView();

            plot.Size = new Size(500, 500);

            PlotModel plotModel = new PlotModel();

            var series = new LineSeries();

            series.Points.AddRange(dataPoints);

            series.LineLegendPosition = LineLegendPosition.Start;

            series.DataFieldX = "X";

            series.DataFieldX = "Y";

            plotModel.Axes.Add(new LogarithmicAxis(){Title = "Error"});
            
            plotModel.Axes.Add(new LinearAxis()
            {
                Title = "Control points",
                Position = AxisPosition.Bottom
                
            });

            plotModel.Series.Add(series);

            plot.Model = plotModel;

            graphForm.Controls.Add(plot);

            graphForm.Size = new Size(700, 700);

            Application.Run(graphForm);
        }
        public static void ShowGraph(Func<double, double> funcSpline, double startx, double endx, double step, string title, List<DataPoint> controlPoints = null, Func<double, double> realFunc = null)
        {

            Form graphForm = new Form();

            PlotView plot = new PlotView();

            plot.Size = new Size(800, 500);

            PlotModel plotModel = new PlotModel();

            plotModel.Series.Add(new FunctionSeries(funcSpline, startx, endx, step, title) {Color = OxyColors.Yellow});

            if (realFunc != null)
            {
                var item = new FunctionSeries(realFunc, startx, endx, step, "Real function") {Color = OxyColors.Purple};
                plotModel.Series.Add(item);
            }


            if (controlPoints != null)
            {
                var pointsSeries = new LineSeries();

                //pointsSeries.LineStyle = LineStyle.None;

                pointsSeries.Points.AddRange(controlPoints);

                pointsSeries.Title = "Control points";

                pointsSeries.Color = OxyColors.Black;

                pointsSeries.MarkerFill = OxyColors.Red;

                pointsSeries.MarkerSize = 5;

                pointsSeries.MarkerType = MarkerType.Circle;

                pointsSeries.Color = OxyColors.Transparent;

                plotModel.Series.Add(pointsSeries);
            }

            plot.Model = plotModel;

            graphForm.Controls.Add(plot);

            graphForm.Size = new Size(900, 700);

            Application.Run(graphForm);
        }
    }
}