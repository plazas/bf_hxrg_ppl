PPL data 

centroid_threshold=0.2
SPOTS: files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-02[2-9][0-9]*.fits")  # andres-02[1-9][0-9]*.fits, 100 ramps of spots, wlamp=180 SPOTS 
FLATS: files1=glob.glob("/projector/aplazas/data/WFIRST/2017-03-02/raw/andres-00[3-9][0-9]_000*.fits")   ## Pick flats away from first ramps to avoid "burn-in", flats at 120W

 flag=False
    if  (np.max(np.abs(res_f))) > 3 or (np.max(np.abs(res_s))) > 3:
        flag=True
        return flag

