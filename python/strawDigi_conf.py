import shipunit as u

class StrawtubesMisalign:
    # TurnOn the module of misalign or no
    misalign = False
    # rand type : "None", "Gaus", "Unif", may added more, need modify strawtubeDigi
    # "None" = all tubes have same maximum sagging, no distribution of max. sag
    # "Gaus" = the tubes have different max sagging with gaus distribution
    # "Unif" = the different max sagging is uniformly distributed in a range given range, not the same as "None"
    randType = "None"

    # maximum sagging at the middle of the tube
    # for using distribution, these are mean value /mpv /something like that
    maxTubeSagging = 0.0*u.cm
    maxWireSagging = 0.0*u.cm

    # uniform distribtion, the delta is half of the range
    tubeUnifDelta = 0.0*u.cm
    wireUnifDelta = 0.0*u.cm

    # Gaussian distribtion, the sigma
    tubeGausSigma = 0.*u.cm
    wireGausSigma = 0.*u.cm

    # debug or not
    debug = False

class DriftTimeCalculate:
    defaultDriftTime = True


