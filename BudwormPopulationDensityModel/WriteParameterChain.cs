using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BudwormPopulationDensityModel
{
    static class WriteParameterChain
    {
        static string  filename;
        public static void AddParameters(int iteration, CalibrationLibrary.Parameters parameters, double period, double p)
        {
            string line = iteration + "\t";
            foreach (CalibrationLibrary.Parameter par in parameters.ModelParameters)
            {
                line += Math.Round(par.RunningValue, 4) + "\t";
            }
            line +=  period + "\t" + Math.Round(p, 4);
            
            System.IO.StreamWriter SW = new System.IO.StreamWriter(filename, true);

            SW.WriteLine(line);

            SW.Close();
        }

        public static void Initialize(string FileName, CalibrationLibrary.Parameters parameters)
        {
            filename = FileName;

            if (System.IO.File.Exists(filename)) System.IO.File.Delete(filename);

            string hdr = "iteration\t";
            foreach (CalibrationLibrary.Parameter p in parameters.ModelParameters)hdr += p.Label + "\t";
            hdr += "Period\tProbability";

            string[] Content = new string[] { hdr };
            System.IO.File.WriteAllLines(FileName, Content);
        }
    }
}
