public static void BlackSholesGetCallValue(double stockPrice, double strikePrice, double riskFreeInterestRate, double timeToMaturity,
                double stdDeviationOfStock, ref double CALL, ref double PUT)
{
    double S = stockPrice; // Stock price
    double K = strikePrice; // Strike price
    double r = riskFreeInterestRate; // risk free interest rate
    double T = timeToMaturity; // time to maturity, days / 365, so this is 3 months
    double stdev = stdDeviationOfStock; // standard deviation of the stock

    double d1 = (Math.Log(S / K) + (r + Math.Pow(stdev, 2) / 2) * T) / (stdev * Math.Sqrt(T));
    double d2 = d1 - stdev * Math.Sqrt(T);

    CALL = S * GetCumulativeNormalDistribution(d1) - K * Math.Exp(-r * T) * GetCumulativeNormalDistribution(d2);
    PUT = K * Math.Exp(-r * T) * GetCumulativeNormalDistribution(-d2) - S * GetCumulativeNormalDistribution(-d1);
}

public static double GetCumulativeNormalDistribution(double x)
{
    double L, K, w;
    /* constants */
    const double a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
    const double a4 = -1.821255978, a5 = 1.330274429;

    L = Math.Abs(x);
    K = 1.0 / (1.0 + 0.2316419 * L);
    w = 1.0 - 1.0 / Math.Sqrt(2 * Math.PI) * Math.Exp(-L * L / 2) * (a1 * K + a2 * K * K + a3 * Math.Pow(K, 3) + a4 * Math.Pow(K, 4) + a5 * Math.Pow(K, 5));

    if (x < 0)
    {
        w = 1.0 - w;
    }
    return w;
}