using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MarkovChain
{
    public class MCMC
    {
        static Random R = new Random();
        static bool accepted;
        
        public static double fractionofdomain;

        public delegate double GetProbability(Parameters parameters);

        
        public static void RunMCMC(Parameters parameters, GetProbability getmodelprobability, int Iterations, double JumpSize)
        {
            fractionofdomain = JumpSize;

            double LogP_last = double.MinValue;

            for (int i = 0; i < Iterations; i++)
            {
                double logp = getmodelprobability(parameters);

                double LogMetropolisAlpha = logp - LogP_last;

                accepted = LogMetropolisAlpha > Math.Log(R.NextDouble());

                if (accepted)
                {
                    parameters.AcceptRunningValues();
                    parameters.Jump(fractionofdomain);
                }
                else
                {
                    parameters.UseLastAcceptedValues();
                }

                LogP_last = logp;
            }
            System.Console.WriteLine("ready");
            //System.Console.ReadLine();
        }

    }
}
