

def spread_pixels(Nside_low, Nside_high, ID):

    from math import log
    Llow = int(log(Nside_low, 2))
    Lhigh = int(log(Nside_high, 2))

    print(Llow, Lhigh)

    b = bin(ID)

    DN = Lhigh-Llow
    a = [bin(i)[2:].zfill(2**DN) for i in range(4**DN)]
    pix_IDs = []
    for i in a:
        x = (b[2:].zfill(Llow) + i)
        pix_IDs.append(int(x, 2))
    
    return(pix_IDs)
