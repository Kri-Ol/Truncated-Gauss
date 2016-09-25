import math
import random
import scipy.special
import scipy.optimize

PI    = 3.1415926535897932384626433832795
SQRT2 = 1.4142135623730950488016887242097

def Phi(z):
    """
    The cumulative function of the unit gaussian distribution
    """

    return  0.5 * (1.0 + scipy.special.erf(z / SQRT2))

def phi(z):
    """
    PDF of the standard gaussian distribution
    """

    return (1.0 / math.sqrt(2.0*PI)) * math.exp(-0.5 * z * z)

def SampleTG(mu, sigma, a, b):
    """
    Sample from
    """

    if sigma <= 0.0:
        raise ValueError("Sample: non-positive sigma")

    # for a moment, work with true bounded gaussian
    if math.isinf(a):
        raise ValueError("Sample: infinite a")
    if math.isinf(b):
        raise ValueError("Sample: infinite b")

    if a >= b:
        raise ValueError("Sample: a >= b")

    while True:
        r = random.gauss(mu, sigma)
        if r >= a and r <= b:
            return r

    return math.Inf

def Mean(mu, sigma, a, b):
    """
    Mean for truncated gaussian
    """

    if sigma <= 0.0:
        raise ValueError("Mean: non-positive sigma")

    # for a moment, work with true bounded gaussian
    if math.isinf(a):
        raise ValueError("Mean: infinite a")
    if math.isinf(b):
        raise ValueError("Mean: infinite b")

    if a >= b:
        raise ValueError("Mean: a >= b")

    alfa = (a - mu) / sigma
    beta = (b - mu) / sigma

    Z = Phi(beta) - Phi(alfa)

    return mu + sigma*(phi(alfa) - phi(beta))/Z

def f(mu, mean, sigma, a, b):
    """
    Function to search for true mu when particular mean is requested
    Root of this function would be right mu
    """
    return mean - Mean(mu, sigma, a, b)

if __name__ == "__main__":

    random.seed(12345)

    # some test printouts
    # print(phi(0.0))
    # print(Phi(0.0))

    a = 50000.0
    b = 250000.0
    mean = 70000.0
    sigma = 24000.0

    mu = scipy.optimize.brentq(f, a, b, args=(mean, sigma, a, b))
    print("Found mu = {0} for the desired mean {1} and sigma {2}".format(mu, mean, sigma))

    # test sampling

    N  = 100000
    s  = 0.0
    s2 = 0.0
    for k in range(0, N):
        q   = SampleTG(mu, sigma, a, b)
        if q < a:
            raise ValueError("Test: sampled value below a")
        if q > b:
            raise ValueError("Test: sampled value above b")
        s  += q
        s2 += q*q

    print("Sampled {0} truncated gaussians and got observed mean = {1}".format(N, s/float(N)))
