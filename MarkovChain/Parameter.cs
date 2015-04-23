using System;
using System.Collections.Generic;
using System.Text;


namespace MarkovChain
{
    
    public class Parameter
    {
        static Random random = new Random();
        //double ShaveFraction;
        public UniversalDistribution distribution;
        List<double> values;
        List<double> probabilities;
        double runningvalue;
        double LastAcceptedValue;
        string label;
        public string Label
        {
            get
            {
                return label;
            }
        }
        public static void SetRandomSeed(int seed)
        {
            random = new Random(seed);
        }
        public double RunningValue
        {
            get
            {
                return runningvalue;
            }
            set
            {
                runningvalue = value;
            }
        }
        double Min(List<double> v)
        {
            double min = double.MaxValue;
            foreach (double _v in v) if (_v <min ) min = _v;
            return min;
        }
        double Max(List<double> v)
        {
            double max = double.MinValue;
            foreach(double _v in v)if(_v > max)max = _v;
            return max;
        }
        
        public void UseLastAcceptedValue()
        {
            runningvalue = LastAcceptedValue;
        }
        public void AcceptRunningValue()
        {
            LastAcceptedValue = runningvalue;
        }
       
        public bool OutOfRange()
        {
            bool ToLarge = (runningvalue > distribution.Max ? true : false);
            bool ToSmall = (runningvalue < distribution.Min ? true : false);
            return (ToLarge || ToSmall ? true : false);
        }
        
        public void Jump(double FractionOfDomain)
        {
            runningvalue = LastAcceptedValue + (random.NextDouble() - 0.5F) * FractionOfDomain * distribution.Range;
        }
        public void RandomJump()
        {
            runningvalue = distribution.Min + random.NextDouble() * distribution.Range;
        }
        public void AddValueAndProbability(double probability)
        {
            values.Add(runningvalue);
            probabilities.Add(probability);

           
        }
        public Parameter(Parameter p)
        {
            this.label = p.label;
            distribution = new UniversalDistribution(p.distribution.Min, p.distribution.Max);
            values = new List<double>();
            probabilities = new List<double>();
        }
        public Parameter(string label, double Min, double Max, double ShaveFraction)
        {
            this.label = label;
            distribution = new UniversalDistribution(Min, Max);
            
            values = new List<double>();
            probabilities = new List<double>();

            // Take a random value within the domain
            RandomJump();
            LastAcceptedValue = runningvalue;
        }
    }
}
