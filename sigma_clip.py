from sys import stdout, stderr
try:
    import numpy
    from numpy import zeros, sqrt
    have_numpy = True
except:
    have_numpy = False


def sigma_clip(arrin, niter=4, nsig=4, get_indices=False, extra={},
               verbose=False, silent=False):
    """
    NAME:
      sigma_clip()

    PURPOSE:
      Calculate the mean/stdev of an array with sigma clipping. Iterate
      niter times, removing elements that are outside nsig, and recalculating
      mean/stdev.

    CALLING SEQUENCE:
      mean,stdev = sigma_clip(arr, niter=4, nsig=4, extra={})

    INPUTS:
      arr: A numpy array or a sequence that can be converted.

    OPTIONAL INPUTS:
      niter: number of iterations, defaults to 4
      nsig: number of sigma, defaults to 4
      get_indices: bool,optional
        if True return mean,stdev,indices

    OUTPUTS:
      mean,stdev: A tuple containing mean and standard deviation.
    OPTIONAL OUTPUTS
      extra={}: Dictionary containing the array of used indices in
         extra['index']

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU
      Minor bug fix to error messaging: 2010-05-28. Brian Gerke, SLAC
      Added silent keyword, to shut off error messages.  BFG 2010-09-13

    """
    arr = numpy.array(arrin, ndmin=1, copy=False)

    index = numpy.arange(arr.size)

    if get_indices:
        res = [None, None, None]
    else:
        res = [None, None]

    for i in numpy.arange(niter):
        m = arr[index].mean()
        s = arr[index].std()

        if verbose:
            stdout.write('iter %s\tnuse: %s\tmean %s\tstdev %s\n' %
                         (i+1, index.size, m, s))

        clip = nsig*s

        w, = numpy.where((numpy.abs(arr[index] - m)) < clip)

        if (w.size == 0):
            if (not silent):
                stderr.write("nsig too small. Everything clipped on "
                             "iteration %d\n" % (i+1))
            res[0] = m
            res[1] = s
            return res

        index = index[w]

    # Calculate final stats
    amean = arr[index].mean()
    asig = arr[index].std()

    res[0] = m
    res[1] = s
    extra['index'] = index
    if get_indices:
        res[2] = index

    return res


if __name__ == "__main__":
    main()
