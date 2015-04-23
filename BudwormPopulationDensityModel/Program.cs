using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
 

namespace BudwormPopulationDensityModel
{
    class Program
    {
        static Random R = new Random();
        static int iteration = 0;
        
        static List<double> Meas;
        
        static double amplitude_opt;

       
        static double LogPmax = double.MinValue;

        public static double GetAbsMax(List<double> Meas)
        {
            double max = double.MinValue;
            foreach (double d in Meas) if (Math.Abs(d) > max) max = d;
            return max;
        }
        public static double RandomGaussian(double mean, double stdDev)
        {
            double u1 = R.NextDouble(); //these are uniform(0,1) random doubles
            double u2 = R.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *Math.Sin(2.0 * Math.PI * u2);  
            double randNormal = mean + stdDev * randStdNormal; //random normal(mean,stdDev^2)

            return randNormal;
        }

        static void PrintToFile(List<double> PopulationSize, string FileName)
        {
            List<string> Content = new List<string>();
            for (int t = 0; t < PopulationSize.Count; t++)
            {
                Content.Add(t + "\t" + PopulationSize[t]);
            }

            System.IO.File.WriteAllLines(FileName, Content.ToArray());
        }
        static double[] RunModel(int NrOfMeas, CalibrationLibrary.Parameters parameters, ref List<double> Model)
        {
            
            double a1 = parameters["a1"].RunningValue;
            double a2 = parameters["a2"].RunningValue;
         
            // array of population size to populate
            Model = new List<double>(new double[NrOfMeas+2]);

            // irruption period
            double d = 1 + a1;
            double e = 4 * a2;

            double w = Math.Atan(Math.Sqrt((Math.Abs(Math.Pow(d, 2) + e))) / d);   //#Royama Eq. 2.17
            double period = 2 * Math.PI / w;                                        //#Royama Eq. 2.19


            // initalize first two years
            Model[0] = 0.4987747;// RandomGaussian(0, 1);
            Model[1] = 0.6088884;// RandomGaussian(0, 1);

            // cycle through years and get simulated pop. size
            double amplitude = 0;
            for (int t = 2; t < NrOfMeas+2; t++)
            {

                Model[t] = (1 + a1) * Model[t - 1] + a2 * Model[t - 2];//  #Royama Eq. 2.6b

                if (Model[t] > amplitude) amplitude = Model[t];
            }

            

            return new double[] { period, amplitude };
        }
        
        static List<double> ReadMeasurements(string FileName)
        {
            
            List<double> Meas = new List<double>();

            string[] Content = System.IO.File.ReadAllLines(FileName);

            for (int line = 1; line < Content.Length;line++ )
            {
                string[] terms = Content[line].Split('\t');
                Meas.Add(double.Parse(terms[1]));
            }

            return Meas;
        }
        static void PrintModelMeas(List<double> Model,List<double> Meas, string FileName)
        {
            List<string> Content = new List<string>();
            Content.Add("Time\tModel\tMeas");

            for (int t = 0; t < Meas.Count; t++)
            {
                Content.Add(t + "\t" + Model[t] + "\t" + Meas[t]);
            }
            System.IO.File.WriteAllLines(FileName, Content.ToArray());
        }
        static double GetModelMeasProb(List<double> Model, List<double> Meas)
        {
            double RSS = 0;
            for (int d = 0; d < Model.Count; d++)
            {
                RSS += Math.Pow(Model[d] - Meas[d], 2);
            }
            return 1.0 / RSS;
        }

        static double GetProbability( CalibrationLibrary.Parameters parameters)
        {

            List<double> Model = null;

            double[] result = RunModel(Meas.Count, parameters, ref Model);


            PrintModelMeas(Model, Meas, "output/compare" + iteration.ToString() + ".txt");

            //double p1 = 1.0 / Math.Pow(result[0] - 5.53, 2);
            double p2 = 1.0/ Math.Pow(result[1] - amplitude_opt, 2);
            double p3 = GetModelMeasProb(Model, Meas);

            double logp =   Math.Log(p3);// +Math.Log(p2); Math.Log(p2) 

            double a1 = Math.Round ( parameters["a1"].RunningValue,2);
            double a2 =  Math.Round (parameters["a2"].RunningValue,2);
        
            double period =  Math.Round ( result[0],2);

            if(iteration >1000)WriteParameterChain.AddParameters(iteration, parameters, period, logp);
            
           // System.Console.WriteLine(period);
            if (logp > LogPmax)
            {
                //System.Console.WriteLine(iteration + "\tp\t" + p + " conv\t" + parameters.GetConvergence());
                string Line = parameters.Values(2);

                System.Console.WriteLine(iteration +"\t" + Math.Round(result[0], 2) + "\t" + Math.Round(result[1], 2) + "\t" +  Math.Round(logp,2));

                System.IO.File.WriteAllLines("output\\BestPV.txt", new string[]{Line});

                PrintModelMeas(Model, Meas, "output/compare.txt");
                LogPmax = logp;
            }
            
            iteration++;

            return logp;
        }
        static bool CheckValidParameterVector(CalibrationLibrary.Parameters parameters)
        {
            foreach (CalibrationLibrary.Parameter p in parameters.ModelParameters)
            {
                if (p.OutOfRange()) return false;
            }

            double a1 = parameters["a1"].RunningValue;
            double a2 = parameters["a2"].RunningValue;

            // irruption period
            double d = 1 + a1;
            double e = 4 * a2;

            double w = Math.Atan(Math.Sqrt((Math.Abs(Math.Pow(d, 2) + e))) / d);   //#Royama Eq. 2.17
            double period = 2 * Math.PI / w;                                        //#Royama Eq. 2.19

            bool Success = (period > 0 ? true : false);

            return Success;

        }
      
        static void Main(string[] args)
        {
             Meas = ReadMeasurements("input\\TimeSeries.txt");
             amplitude_opt = GetAbsMax(Meas);

            //X = stochasticity (noise.scale)
            //Y = density dependence parameter (a1)
            //Z = delayed density dependence parameter (a2)

            CalibrationLibrary.Parameters parameters = new CalibrationLibrary.Parameters(CheckValidParameterVector);

            //parameters.Add(new CalibrationLibrary.Parameter("noise_scale",1,1, 0.5));
            //parameters.Add(new CalibrationLibrary.Parameter("a1", -0.9, -0.9, 0.5));
            //parameters.Add(new CalibrationLibrary.Parameter("a2", 0.6, 0.6, 0.5));

            parameters.Add(new CalibrationLibrary.Parameter("a1", -5, 5, 0.5));
            parameters.Add(new CalibrationLibrary.Parameter("a2", -5, 5, 0.5));

            parameters["a1"].RunningValue = 0.4; 
            parameters["a2"].RunningValue = 0.3; 
           

            WriteParameterChain.Initialize("output/ParameterEstimation.txt", parameters);

            //CalibrationLibrary.RunIterations r = new CalibrationLibrary.RunIterations(parameters, int.MaxValue, GetProbability);

            //CalibrationLibrary.Calibration c = new CalibrationLibrary.Calibration(parameters, GetProbability, int.MaxValue, 1000, 0.1);

            CalibrationLibrary.MCMC.RunMCMC(parameters, GetProbability, 500000, 0.1F);

           
       
        }
    }
}
