def calculate(aaseq):
    hydrolevels=[]
    windowscores=[]
    hydroscore=0
    counter=0
    lowest=7000
    hydrophob = {
            'S':    0.84,   # Serine
            
            'F' : -0.32,    # Phenylalanine
            
            'L' : -0.55,    # Leucine
             
            'Y' : 0.68,    # Tyrosine
             
            '_' : 0,    # Stop
             
            'C' : -.13,    # Cysteine
             
            'W' : 0.30,    # Tryptophan
            
             
            'P' : 2.23,    # Proline
            
            'H' : 2.06,    # Histidine
            
            'Q' : 2.36,    # Glutamine
            
            'R' : 2.58,    # Arginine
            
            'I' : -0.60,    # Isoleucine
            
            'M' : -0.10,    # Methionine
            'T' : 0.52,    # Threonine
            'N' : 2.05,    # Asparagine
            'K' : 2.71,    # Lysine
            'S' : 0.84,    # Serine
            'R' : 2.58,    # Arginine
          
            'V' : -0.31,    # Valine
          
            'A' : 0.11,    # Alanine
          
            'D' : 3.49,    # Aspartic Acid
            'E' : 2.68,    # Glutamic Acid
            'G' : 0.74,    # Glycine
                }

    for i in range(0, len(aaseq)-10):
        hydroscore= 0
        for w in range(10):
            hydroscore += hydrophob[aaseq[i+w]]
        if hydroscore<lowest :
            lowest = hydroscore
        windowscores.append(hydroscore)
        ##hydrolevels.append(hydrophob[aaseq[i]])


        if hydrophob[aaseq[i]] == '_':
            break

    # Return
    print(lowest)
    return 0
    #return windowscores
