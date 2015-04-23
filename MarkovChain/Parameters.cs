using System;
using System.Collections.Generic;
using System.Text;

namespace MarkovChain
{
    public delegate bool ValidParameterVector(Parameters parameters);

    public class Parameters
    {
        ValidParameterVector checkparameters;

        List<Parameter> modelparameters;

        public bool checkoutofrange(Parameters parameters)
        {
            foreach (Parameter p in parameters.ModelParameters)
            {
                if (p.OutOfRange()) return false;
            }
            return true;
        }

         

        public List<Parameter> ModelParameters
        {
            get
            {
                return modelparameters ;
            }
        }

        public Parameter this[string label]
        {
            get
            {
                return Get(label);
            }
        }
        public void Add(Parameter p)
        {
            modelparameters.Add(p);
        }
        public string Values(int decimals)
        {
            string line = "";
            foreach (Parameter p in ModelParameters)
            {
                line += p.Label + "\t" + Math.Round(p.RunningValue,2) + "\t";
            }
            return line;
        }
        private Parameter Get(string label)
        {
            foreach (Parameter p in ModelParameters)
            {
                if (p.Label == label) return p;
            }
            return null;          
        }
        public void UseLastAcceptedValues()
        {
            foreach (Parameter p in ModelParameters)
            {
                p.UseLastAcceptedValue();
            }
        }
        public void AcceptRunningValues()
        {
            foreach (Parameter p in ModelParameters)
            {
                p.AcceptRunningValue();
            }
        }
        
        public void Jump(double FractionOfDomain)
        {
            
            do
            {
                foreach (Parameter p in ModelParameters)
                {
                    p.Jump(FractionOfDomain);
                }
            }
            while (checkparameters(this) == false);
        }
        public void RandomJump()
        {
            do
            {
                foreach (Parameter p in ModelParameters)
                {
                    p.RandomJump();
                }
            }
            while (!checkparameters(this));
        }

        public Parameters(ValidParameterVector Checkparameters)
        {
            checkparameters = Checkparameters;
            modelparameters = new List<Parameter>();
        }
        public Parameters()
        {
            checkparameters = checkoutofrange;
            modelparameters = new List<Parameter>();
        }
    }
}
