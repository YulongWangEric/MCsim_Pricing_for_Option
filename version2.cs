/**
 * This package implements a Project using MC simulation to price European Call/Put Options
 * under Euler Discretization scheme
 * 
 * @author Yulong Wang
 * @version 02/09/15
 */

using System;

namespace OptionPricing
{
    public enum OptionType { Call, Put };

    ///   <summary> 
    ///   This class generates Gaussian distributed pseudo-random variates 
    ///   </summary> 

    public class GaussianGenerator
    {
        private Random rand;

        public void SetSeed(int seed)
        {
            rand = new Random(seed);
        }

        ///   <summary> 
        ///   this method uses two U(0,1) distributed random variates 
        /// to produce a Gaussian distributed random variate
        ///   </summary> 

        public double NextGaussian()
        {
            double u1 = rand.NextDouble();
            double u2 = rand.NextDouble();
            double u = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
            return u;
        }

        ///   <summary> 
        ///   constructor with no parameter. A random integer generated 
        /// by system time is employed as the seed
        ///   </summary> 

        public GaussianGenerator()
        {
            int seed = (new Random()).Next();
            rand = new Random(seed);
        }

        /// <summary>
        /// overloading constructor given an integer parameter
        /// </summary>
        /// <param name="seed"> is used to set the seed for<member>rand</member></param>

        public GaussianGenerator(int seed)
        {
            rand = new Random(seed);
        }
    }


    ///   <summary> 
    ///   This is an abstract class for the object of plain vanilla option 
    ///   </summary> 

    public abstract class PlainVanillaOption
    {
        /// <summary>
        /// <member>t is the time to maturity for this option</member>
        /// <member>strikePrice is the strike price for this option</member>
        /// <member>underlyingAssetName is the name for underlying asset</member>
        /// <member>oType defines whether it is a call or put option</member>
        /// </summary>

        protected double t;
        protected double strikePrice;
        protected string underlyingAssetName;
        protected OptionType oType;

        public double T
        {
            get { return t; }
            set { t = value; }
        }

        public double StrikePrice
        {
            get { return strikePrice; }
            set { strikePrice = value; }
        }

        public string UnderlyingAssetName
        {
            get { return underlyingAssetName != null ? underlyingAssetName : "NA"; }
            set { underlyingAssetName = value; }
        }

        public OptionType Otype
        {
            get { return oType; }
        }

        /// <summary>
        /// method which returns the value of an option at maturity 
        /// given its underlying asset price at maturity
        /// </summary>
        /// <param name="finalPrice"> assumes the underlying asset price at maturity</param>
        /// <returns> A double type value which describes the future value of an option </returns>

        public double ValueAtMaturity(double finalPrice)
        {
            if (oType == OptionType.Call)
            { return Math.Max(0.0, finalPrice - strikePrice); }
            else { return Math.Max(0.0, strikePrice - finalPrice); }
        }

        /// <summary>
        /// This method gives the estimated the option value 
        /// and the estimation standard error using MC simulation
        /// </summary>
        /// <param name="S"> describes the Brownian motion rule of the underlying asset</param>
        /// <param name="D"> chooses the discretization scheme </param>
        /// <param name="numberOfScenarios">decides the number of random paths to be generated</param>
        /// <param name="timeSteps"> decides the number of intervals for discretization scheme </param>
        /// <param name="antitheticFlag"> chooses whether to use antithetic variance reduction techniques</param>
        /// <returns> a 2 elements array. The first one gives the estimated option price, 
        /// the second gives the estimation error</returns>

        public abstract double[] PricingByMCSim(StochasticAssetPrice S, IDiscretizationScheme D,
            int numberOfScenarios, int timeSteps, bool antitheticFlag);
    }


    /// <summary>
    /// This class describes the Brownian Motion rule for the price of an asset
    /// </summary>

    public class StochasticAssetPrice
    {
        /// <summary>
        /// <member> mu is the drift parameter</member>
        /// <member> sigma is the volatility parameter</member>
        /// <member> currentPrice is the underlying asset spot price</member>
        /// </summary>

        private double mu;
        private double sigma;
        private double currentPrice;

        public double Mu
        {
            get { return mu; }
            set { mu = value; }
        }

        public double Sigma
        {
            get { return sigma; }
            set { sigma = value; }
        }

        public double CurrentPrice
        {
            get { return currentPrice; }
            set { currentPrice = value; }
        }

        /// Constructor
        public StochasticAssetPrice(double Mu, double sigma, double currentPrice)
        {
            this.mu = Mu;
            this.sigma = sigma;
            this.currentPrice = currentPrice;
        }

        /// Copy method
        public StochasticAssetPrice(StochasticAssetPrice S)
        {
            this.Mu = S.Mu;
            this.sigma = S.Sigma;
            this.currentPrice = S.CurrentPrice;
        }
    }


    /// <summary>
    /// This interface defines the framework for discretization scheme
    /// </summary>

    public interface IDiscretizationScheme
    {
        /// <summary>
        /// This method simulates a one step motion for the underlying
        /// asset price based on specified discretization scheme and time interval 
        /// </summary>
        /// <param name="S">reference to the current state of underlying asset </param>
        /// <param name="dt"> decides the time length for each discretized interval</param>
        /// <param name="Z"> is the random factors for Brownian motion</param>
        /// <returns> the asset price at the next step </returns>

        double GetNextPrice(StochasticAssetPrice S, double dt, params double[] Z);

        /// <summary>
        /// This method simulates a whole path for the underlying asset 
        /// price based on specified discretization scheme and time to maturity 
        /// </summary>
        /// <param name="S">reference to the current state of underlying asset </param>
        /// <param name="nrandom">reference to the Gaussian random variates generator</param>
        /// <param name="totalTime">defines the time to maturity</param>
        /// <param name="timeSteps">decides the number of discretized intervals for the whole period</param>
        /// <returns>an double type array which records the asset price at each point</returns>

        double[] GeneratingRandomPricePath(StochasticAssetPrice S, GaussianGenerator nrandom,
            double totalTime, int timeSteps);

        /// <summary>
        /// This method simulates two price path using two series of mutual negative 
        /// random variates named by 'antithetic' which can help reduce the variance 
        /// of MC simulation
        /// </summary>
        /// <returns>an double type 2*N array which records the asset price at each point
        /// for the two paths
        /// </returns>

        double[][] DipathByAntitheticMethod(StochasticAssetPrice S, GaussianGenerator nrandom,
            double totalTime, int timeSteps);
    }

    /// <summary>
    /// This abstract class implements the IDiscretizationScheme interface
    /// based on Black Scholes Model (mu and sigma is constant) 
    /// </summary>

    public abstract class DiscretizationSchemeForBSModel : IDiscretizationScheme
    {

        public abstract double GetNextPrice(StochasticAssetPrice S, double dt, params double[] Z);

        /// <summary>
        /// override the <method>GeneratingRandomPricePath</method>.
        /// For Black Scholes model, mu and sigma is constant, and only one random
        /// variate is needed to generate the state for next point.
        /// </summary>

        public double[] GeneratingRandomPricePath(StochasticAssetPrice S, GaussianGenerator nrandom,
            double totalTime, int timeSteps)
        {
            double[] pricePath = new double[timeSteps + 1];
            double dt = totalTime / (double)timeSteps;
            pricePath[0] = S.CurrentPrice;
            for (int i = 1; i <= timeSteps; i++)
            {
                double z = nrandom.NextGaussian();
                pricePath[i] = GetNextPrice(S, dt, z);
            }
            return pricePath;
        }

        /// <summary>
        /// override the <method>DipathByAntitheticMethod</method>.
        /// For Black Scholes model, mu and sigma is constant, and only one random
        /// variate is needed to generate the state for next point
        /// </summary>

        public double[][] DipathByAntitheticMethod(StochasticAssetPrice S, GaussianGenerator nrandom,
            double totalTime, int timeSteps)
        {
            double[][] pricePath = new double[2][];
            pricePath[0] = new double[timeSteps + 1];
            pricePath[1] = new double[timeSteps + 1];
            StochasticAssetPrice S2 = new StochasticAssetPrice(S);
            double dt = totalTime / (double)timeSteps;
            pricePath[0][0] = S.CurrentPrice;
            pricePath[1][0] = S.CurrentPrice;
            for (int i = 1; i <= timeSteps; i++)
            {
                double z = nrandom.NextGaussian();
                double z2 = -z;
                pricePath[0][i] = GetNextPrice(S, dt, z);
                pricePath[1][i] = GetNextPrice(S2, dt, z2);
            }
            return pricePath;
        }
    }

    /// <summary>
    /// This class inherites the DiscretizationSchemeForBSModel abstract class
    /// by specify the calculation method for the next stage asset price
    /// </summary>

    public class EulerSchemeForBSModel : DiscretizationSchemeForBSModel
    {
        /// <summary>
        /// Override the <method>GetNextPrice</method> based on Euler Scheme
        /// The StochasticAssetPrice object is updated in this method
        /// </summary>
        /// <returns> a double type value for the price of next stage</returns>

        public override double GetNextPrice(StochasticAssetPrice S, double dt, params double[] Z)
        {
            double nextPrice = S.CurrentPrice * Math.Exp((S.Mu - 0.5 * S.Sigma * S.Sigma) * dt
                + S.Sigma * Z[0] * Math.Sqrt(dt));
            S.CurrentPrice = nextPrice;
            return nextPrice;
        }
    }

    /// <summary>
    /// This class inherites PlainVanillaOption class
    /// The European Option can only be executed at maturity, thus its value does not
    /// depend on the path before maturity
    /// </summary>

    public class EuropeanOption : PlainVanillaOption
    {

        public EuropeanOption(double t, double strikePrice, OptionType oType, string underlyingAssetName = "")
        {
            this.t = t;
            this.oType = oType;
            this.strikePrice = strikePrice;
            this.underlyingAssetName = underlyingAssetName;
        }

        /// <summary>
        /// override the <method>PricingByMCSim</method>
        /// In each scenario, a price path is generated and the option value under 
        /// this scenario is only decided by the final price. So, we take the average 
        /// for the 'value at maturity' of all the simulated scenarios and 
        /// discounted as the approximation for the option value. 
        /// This method also outputs the standard error for the approximation.
        /// </summary>

        public override double[] PricingByMCSim(StochasticAssetPrice S, IDiscretizationScheme D, int numOfScenarios,
            int timeSteps, bool antitheticFlag)
        {
            double[] Result = new double[2];                     //record the approximated price and standard error
            double[] Sample = new double[numOfScenarios];        //record value at maturity for each scenario
            double sumValue = 0.0;                               //record cumulated sum of Sample array
            GaussianGenerator nrand = new GaussianGenerator();
            if (antitheticFlag)
            {
                //when antithetic variance reduction technique is used
                for (int i = 0; i < numOfScenarios; i++)
                {
                    StochasticAssetPrice S1 = new StochasticAssetPrice(S);
                    double[][] PricePath;
                    PricePath = D.DipathByAntitheticMethod(S1, nrand, T, timeSteps);
                    Sample[i] = (ValueAtMaturity(PricePath[0][timeSteps])
                        + ValueAtMaturity(PricePath[1][timeSteps])) / 2.0;
                    sumValue += Sample[i];
                }
            }
            else
            {
                //when antithetic variance reduction technique is not used
                for (int i = 0; i < numOfScenarios; i++)
                {
                    StochasticAssetPrice S1 = new StochasticAssetPrice(S);
                    double[] PricePath;
                    PricePath = D.GeneratingRandomPricePath(S1, nrand, T, timeSteps);
                    Sample[i] = ValueAtMaturity(PricePath[timeSteps]);
                    sumValue += Sample[i];
                }
            }
            double disFactor = Math.Exp(-S.Mu * T);             //discount factor to convert value at maturity to current
            Result[0] = sumValue / (double)numOfScenarios * disFactor;
            double totalVariance = 0.0;
            for (int i = 0; i < numOfScenarios; i++)
            {
                totalVariance += Math.Pow(Sample[i] - Result[0], 2.0);
            }
            if (numOfScenarios > 1)
            {
                double std = Math.Sqrt(totalVariance / (numOfScenarios - 1));
                Result[1] = std / Math.Sqrt(numOfScenarios) * disFactor;
            }
            return Result;
        }
    }

    /// <summary>
    /// testing program
    /// </summary>

    class Program
    {
        static void Main(string[] args)
        {
            //input needed information from keyboard
            Console.WriteLine("Input spot price of underlying asset:");
            double spotPrice = Convert.ToDouble(Console.ReadLine());
            int oType;
            do
            {
                Console.WriteLine("Input European option type (0 for call, 1 for put):");
                oType = Convert.ToInt32(Console.ReadLine());
            }
            while (oType != 0 & oType != 1);
            OptionType optionType = (oType == 0 ? OptionType.Call : OptionType.Put);
            Console.WriteLine("Input option strike price:");
            double strikePrice = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Input time to maturity of this option:");
            double timeToMaturity = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Input drift parameter for brownian motion:");
            double Mu = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Input volatility parameter for brownian motion:");
            double sigma = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Input number of scenarios generated by MC simulation:");
            int numOfScenarios = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Input number of time steps for Euler discretization:");
            int timeSteps = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Input whether to use antithetic variance reduction technique " +
                "(true for yes,false for no): ");
            bool antithetic = Convert.ToBoolean(Console.ReadLine());
            EuropeanOption Option = new EuropeanOption(timeToMaturity, strikePrice, optionType);
            EulerSchemeForBSModel Euler = new EulerSchemeForBSModel();
            StochasticAssetPrice Asset = new StochasticAssetPrice(Mu, sigma, spotPrice);
            if (antithetic)
            {
                double[] s1 = Option.PricingByMCSim(Asset, Euler, numOfScenarios, timeSteps,
                    true);
                Console.WriteLine("Option price estimated by Monte Carlo Simulation and " +
                    "antithetic variance reduction is:\n {0:#0.00} \n Standard error is:\n {1:#0.000} ",
                    s1[0], s1[1]);
            }
            else
            {
                double[] s2 = Option.PricingByMCSim(Asset, Euler, numOfScenarios, timeSteps, true);
                Console.WriteLine("Option price estimated by Monte Carlo Simulation is:\n {0:#0.00} \n" +
                    "Standard error is:\n {1:#0.000} ", s2[0], s2[1]);
            }
            Console.ReadLine();
        }
    }
}
