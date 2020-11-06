
import matplotlib
matplotlib.use('TKagg')
import matplotlib.pyplot as plt
import numpy as np
import sys

from multi_isotope_calculator import Multi_isotope


def main():
    #compare_alpha_centrifuge()
    #compare_nat_centrifuge()
    #compare_reprocessed_centrifuge()
    
    #compare_alpha_diffusion()
    #compare_nat_diffusion()
    compare_reprocessed_diffusion()

    return

def comparison(wood, multi, string, latex=False):
    """PPrint the parameters. If latex=True, print latex-table style"""
    print_opt = np.get_printoptions()
    np.set_printoptions(precision=5)
    diff = 100 * (wood-multi) / wood
    
    print('\n\nComparison of ' + string)

    if latex:
        for i in range(len(wood)):
            print("&{:.5e}  &{:.5e} &{:.5f} \\\\".format(multi[i], 
                                                         wood[i], 
                                                         diff[i]))
        np.set_printoptions(print_opt)
        return

    print('Wood:                 ', wood)
    print('Own calculations:     ', multi)
    print("Diff [%]:             ", (wood-multi)/wood)
    
    return

def compare_alpha_centrifuge():
    c = Multi_isotope({'235': (0.72, 3, 0.3)}, process='centrifuge', feed=1)
    wood = np.array([2.2, 2.0, 1.8, 1.6, 1.4, 1.])
    multi = np.array(c.alpha)
    
    comparison(wood, multi, 'stage separation factors') 
    return

def compare_nat_centrifuge():
    wood = np.array([[0.045739, 5], [0.889859, 93]])
    
    natural = {'234': 5.5e-3, '235': (0.72, 5, 0.3)}
    c = Multi_isotope(natural, process='centrifuge', product=1)
    
    for i, p in enumerate(wood[:,1]):
        c.set_product_enrichment(p)
        c.calculate_staging()
    
        comparison(wood[i], c.xp[2:4]*100, 
                   ('product U-234 concentration [%]\n'
                    + 'natural U to {} %'.format(p)))
    return

def compare_reprocessed_centrifuge():
    wood = np.array([[7.2635e-9, 0.128628, 5, 1.673063], 
                     [1.1531e-7, 2.016983, 75, 19.342117]])
    
    reprocessed = {'232': 1e-9, '234': 0.02, '235': (0.9, 5, 0.3), 
                    '236': 0.4}
    c = Multi_isotope(reprocessed, process='centrifuge', product=1)
    
    for i, p in enumerate(wood[:,2]):
        c.set_product_enrichment(p)
        c.calculate_staging()
    
        comparison(wood[i], c.xp[[0,2,3,4]]*100, 
                   ('product U232, 234, 235, 236 concentrations [%]\n'
                    + 'reprocessed U to {} %'.format(p)), latex=True)
    return

def compare_alpha_diffusion():
    c = Multi_isotope({'235': (0.72, 3, 0.3)}, process='diffusion', feed=1)
    wood = np.array([1.008633, 1.007179, 1.005731, 1.004289, 1.002853, 1])
    multi = np.array(c.alpha)
    
    comparison(wood, multi, 'stage separation factors', latex=True)
    return

def compare_nat_diffusion():
    wood = np.array([[0.047791, 5], [0.934694, 93]])
    
    natural = {'234': 5.5e-3, '235': (0.72, 5, 0.3)}
    c = Multi_isotope(natural, process='diffusion', product=1)
    
    for i, p in enumerate(wood[:,1]):
        c.set_product_enrichment(p)
        c.calculate_staging()
    
        comparison(wood[i], c.xp[2:4]*100, 
                   ('product U-234 concentration [%]\n'
                    + 'natural U to {} %'.format(p)), latex=True)
    return

def compare_reprocessed_diffusion():
    wood = np.array([[7.5567e-9, 0.133083, 5, 1.564999], 
                     [1.2008e-7, 2.097610, 75., 16.348610]])
    
    reprocessed = {'232': 1e-9, '234': 0.02, '235': (0.9, 5, 0.3), 
                    '236': 0.4}
    c = Multi_isotope(reprocessed, process='diffusion', product=1)
    
    for i, p in enumerate(wood[:,2]):
        c.set_product_enrichment(p)
        c.calculate_staging()
    
        comparison(wood[i], c.xp[[0,2,3,4]]*100, 
                   ('product U232, 234, 235, 236 concentrations [%]\n'
                    + 'reprocessed U to {} %'.format(p)), latex=True)
    return

if __name__=='__main__':
    main()
