using Accord.Statistics;
using Accord.Statistics.Distributions.Univariate;
using Accord.Statistics.Testing;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading;

namespace windows_cep_calculator
{
    class Program
    {
        static void Main(string[] args)
        {
            Thread.CurrentThread.CurrentCulture = CultureInfo.CreateSpecificCulture("en-GB");

            if (args.Length == 0)
            {
                Console.WriteLine("Drag a csv file with the results to the exe and try again.");
                Console.Read();
                return;
            }

            string filePath = args[0];
            string[] fileContent = File.ReadAllLines(filePath);

            string resultContent = "";


            resultContent += ("=============================================\n");
            resultContent += ("CEP Calculator\n");
            resultContent += ("Author: Wilson de Miranda, Jean Knapp\n");
            resultContent += ("=============================================\n");

            List<(double, double, double, double)> values = new List<(double, double, double, double)>();

            List<double> column1List = new List<double>();
            List<double> column2List = new List<double>();
            List<double> column3List = new List<double>();
            List<double> column4List = new List<double>();

            for (int i = 1; i < fileContent.Length; i++)
            {
                string line = fileContent[i];

                if (line.Length == 0)
                    continue;

                if (line.Contains(";"))
                {
                    // Está em portugues
                    line = line.Replace(",", ".").Replace(";", ",");
                } else
                {
                    // Está em ingles
                }

                string[] cells = line.Split(',');

                if (cells.Length < 2)
                    continue;
                double a = double.Parse(cells[0]);
                double b = double.Parse(cells[1]);
                double c = Math.Sqrt(a * a + b * b);

                double d;
                if (c == 0)
                {
                    d = Math.Log(1e-1);
                } else
                {
                    d = Math.Log(c);
                }               

                values.Add((a,b,c,d));

                column1List.Add(a);
                column2List.Add(b);
                column3List.Add(c);
                column4List.Add(d);
            }

            double[] column1 = column1List.ToArray();
            double[] column2 = column2List.ToArray();
            double[] column3 = column3List.ToArray();
            double[] column4 = column4List.ToArray();

            // Média amostral
            int n = values.Count;

            double mean1 = Measures.Mean(column1);
            double mean2 = Measures.Mean(column2);
            double mean3 = Measures.Mean(column3);
            double mean4 = Measures.Mean(column4);

            // Desvio padrão amostral
            double var1 = Measures.Variance(column1, mean1);
            double var2 = Measures.Variance(column2, mean2);
            double var3 = Measures.Variance(column3, mean3);
            double var4 = Measures.Variance(column4, mean4);

            // Mediana
            double median1 = Measures.Median(column1);
            double median2 = Measures.Median(column2);
            double median3 = Measures.Median(column3);
            double median4 = Measures.Median(column4);

            resultContent += ("\nSample:\n");

            var sb = new System.Text.StringBuilder();
            sb.Append(String.Format("{0,10} {1,10} {2,10} {3,10}\n", "", "xbar", "s^2", "median"));
            sb.Append(String.Format("{0,10} {1,10} {2,10} {3,10}\n", "deflection", Math.Round(mean1,2) + "m", Math.Round(var1) + "m", Math.Round(median1, 2) + "m"));
            sb.Append(String.Format("{0,10} {1,10} {2,10} {3,10}\n", "range", Math.Round(mean2) + "m", Math.Round(var2) + "m", Math.Round(median2, 2) + "m"));
            resultContent += (sb + "\n");

            resultContent += ("\nNormality test: \n");
            bool test1a = ShapiroTest(column1);
            bool test1b = LillieTest(column1);
            bool test1c = AndersonTest(column1);
            bool test2a = ShapiroTest(column2);
            bool test2b = LillieTest(column2);
            bool test2c = AndersonTest(column2);

            sb = new System.Text.StringBuilder();
            sb.Append(String.Format("{0,10} {1,10} {2,1} {3,10} {4,1} {5,10}\n", 
                "", 
                "Shapiro",
                "|",
                "Lillie",
                "|",
                "Anderson"));
            sb.Append(String.Format("{0,10} {1,10} {2,1} {3,10} {4,1} {5,10}\n", 
                "range", 
                (test1a ? "Normal" : "Not normal"),
                "|",
                (test1b ? "Normal" : "Not normal"),
                "|",
                (test1c ? "Normal" : "Not normal")));
            sb.Append(String.Format("{0,10} {1,10} {2,1} {3,10} {4,1} {5,10}\n", 
                "deflection",
                (test2a ? "Normal" : "Not normal"),
                "|",
                (test2b ? "Normal" : "Not normal"),
                "|",
                (test2c ? "Normal" : "Not normal")));
            resultContent += (sb + "\n");

            int column1NormalCounter = (test1a ? 1 : 0) + (test1b ? 1 : 0) + (test1c ? 1 : 0);
            int column2NormalCounter = (test2a ? 1 : 0) + (test2b ? 1 : 0) + (test2c ? 1 : 0);

            double CEP = 0;

            if (column1NormalCounter < 2 || column2NormalCounter < 2)
            {
                resultContent += ("\nResult: Not normal\n");
                // CEP Ethridge
                double kurtosis4 = Kurtosis(column4);

                double Dj = 0;
                double[] Di = new double[values.Count];
                double[] T = new double[values.Count];
                double[] Wi = new double[values.Count];
                double[] Ui = new double[values.Count];

                for (int i = 0; i < values.Count; i++)
                {
                    T[i] = 1 - (0.03 * Math.Pow(kurtosis4 - 3, 3) * (Math.Pow(column4[i] - median4, 2) / var4));
                    Di[i] = Math.Max(T[i], 0.01);

                    Dj = Dj + 1 / Di[i];
                }

                double Utotal = 0;
                for (int i = 0; i < values.Count; i++)
                {
                    Wi[i] = 1 / Di[i] / Dj;
                    Ui[i] = Wi[i] * column4[i];

                    Utotal = Utotal + Ui[i];
                }
                
                CEP = Math.Exp(Utotal);

                resultContent += ("\nCEP Equivalent: " + Math.Round(CEP,2) + "m (Ethridge)\n");
            } else
            {
                // Normal
                resultContent += ("\nResult: Normal\n");
                FDistribution F = new FDistribution(degrees1: n -1, degrees2: n - 1);
                double li = var1 / var2 / F.InverseDistributionFunction(p: 0.975);
                double ls = var1 / var2 / F.InverseDistributionFunction(p: 0.025);

                if (li <= 1 && ls >= 1)
                {
                    resultContent += ("        Variances are equal\n");
                    CEP = 1.1774 * (Math.Sqrt(var1) + Math.Sqrt(var2)) / 2;

                    resultContent += ("\nCEP Equivalent: " + Math.Round(CEP, 2) + "m (Closed form integration)\n");
                    resultContent += ("DEP Equivalent: " + Math.Round(Math.Sqrt(var1) * 0.6745, 2) + "m (Closed form integration)\n");
                    resultContent += ("REP Equivalent: " + Math.Round(Math.Sqrt(var2) * 0.6745, 2) + "m (Closed form integration)\n");

                } else
                {
                    resultContent += ("        Variances are different\n");
                    // CEP Grubbs
                    var sigma2 = Math.Abs(var1 + var2);
                    var sigma4 = sigma2 * sigma2;

                    // Aim point
                    double a = 0;
                    double b = 0;

                    double m = 1 + 1 / sigma2 * ((a - mean1) * (a - mean1) + (b - mean2) * (b - mean2));
                    double v = 2 * ((var1 * var1 + var2 * var2) + 2 * (var1 * (mean1 - a) * (mean1 - a) + var2 * (mean2 - b) * (mean2 - b))) / sigma4;

                    CEP = Math.Sqrt(sigma2 * m * (1 - v / (9 * m * m)) * (1 - v / (9 * m * m)) * (1 - v / (9 * m * m)));

                    resultContent += ("\nCEP Equivalent: " + Math.Round(CEP, 2) + "m (Grubbs)\n");

                    // CEP 
                    double cep2 = 0.562 * Math.Max(Math.Sqrt(var1), Math.Sqrt(var2)) + 0.617 * Math.Min(Math.Sqrt(var1), Math.Sqrt(var2));
                    resultContent += ("\nCEP Equivalent: " + Math.Round(cep2, 2) + "m (Pittman)\n");

                    resultContent += ("\nDEP Equivalent: " + Math.Round(Math.Sqrt(var1) * 0.6745, 2) + "m (Closed form integration)\n");
                    resultContent += ("REP Equivalent: " + Math.Round(Math.Sqrt(var2) * 0.6745, 2) + "m (Closed form integration)\n");

                }

            }

            Console.WriteLine(resultContent);
            File.WriteAllText(filePath.Substring(0, filePath.Length - 4) + "_output.txt", resultContent);

            //Console.WriteLine("CEP Equivalente: " + CEP);
            //System.Diagnostics.Debugger.Break();
            Console.ReadLine();


        }

        public static double Kurtosis(double[] sample)
        {
            double mu = Measures.Mean(sample);
            int n = sample.Length;

            double num = 0;
            for (int i = 0; i < sample.Length; i++)
            {
                num = num + Math.Pow(sample[i] - mu, 4) * n;
            }

            double den = 0;
            for (int i = 0; i < sample.Length; i++)
            {
                den = den + Math.Pow(sample[i] - mu, 2);
            }
            den = Math.Pow(den, 2);

            double k = num / den;

            return k;
        }

        public static bool ShapiroTest(double[] sample)
        {
            var result = new ShapiroWilkTest(sample);
            return (result.PValue >= 0.05);
        }

        public static bool LillieTest(double[] sample)
        {
            var result = new KolmogorovSmirnovTest(sample, NormalDistribution.Standard);
            return (result.PValue >= 0.05);
        }

        public static bool AndersonTest(double[] sample)
        {
            NormalDistribution distribution = NormalDistribution.Estimate(sample);
            var result = new AndersonDarlingTest(sample, distribution);
            return (result.PValue >= 0.05);
        }
    }
}
