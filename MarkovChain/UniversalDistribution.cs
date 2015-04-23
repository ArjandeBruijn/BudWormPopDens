using System;
using System.Collections.Generic;
using System.Text;

namespace MarkovChain
{
    public class UniversalDistribution
    {
        public double Min;
        public double Max;

        public bool InDomain(double v)
        {
            if (v > Min && v < Max) return true;
            return false;
        }
        public double Range
        {
            get
            {
                return Max - Min;
            }
        }
        public UniversalDistribution(double Min, double Max)
        {
            this.Min = Min;
            this.Max = Max;
        }
    }
}
