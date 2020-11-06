
import logging
import numpy as np
from scipy.optimize import minimize as scipy_minimize
from sys import stderr as sys_stderr
import warnings

"""
To do list:
    - implement warnings using warnings module

    - check all of the error/exception messages
    
    - add papers to concentration_diffusion and concentration_centrifuge
"""

#logging.basicConfig(stream=sys_stderr, level=logging.DEBUG)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - CLASSES - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class Multi_isotope:
    def __init__(self, concentrations, max_swu=np.inf, n_e=1000, n_s=1000, 
                 process='centrifuge', feed=np.inf, product=np.inf,
                 alpha235=1.6, mstar=350.5, downblend=True):
        """
        Calculate enrichment, cascade parameters for multicomponent isotope
        enrichment.

        This class allows to calculate the enrichment parameters of a 
        cascade of given shape taking into account U-232, U-233, U-234, 
        U-235 and U-238.

       
        References:
        [1] Houston G. Wood, 'Effects of Separation Processes on Minor 
            Uranium Isotopes in Enrichment Cascades'. Science and Global 
            Security, 16:26–36 (2008), DOI: 10.1080/08929880802361796.
        [2] E. von Halle, 'Multicomponent Isotope Separation in Matched 
            Abundance Ratio Cascades Composed of Stages With Large 
            Separation Factors'. International Technology Programs 
            (K/ITP--131). Oak Ridge, TN (1987).
        [3] Richard D. Harvey, 'An Optimization Method for Matched 
            Abundance-Ratio Cascades by Varying the Key Weight'. PhD 
            dissertation. University of Tennessee, Knoxville, TN (2017).
        """
        # If the values (most notably the number of stages and the 
        # enrichment concentrations) are up-to-date, set variable to true.
        # Else, e.g., after changing desired U235 product concentration, 
        # set it to False.
        self.uptodate = False
 
        self.check_input(concentrations, process, alpha235, feed, product,
                         max_swu)
        
        self.mstar = mstar
       
        m_uf6 = 6*19
        self.isotope = (np.array([232, 233, 234, 235, 236, 238], dtype=int)
                        + m_uf6)
        self.xf, self.xp, self.xt = self.init_concentrations(concentrations)
        self.user_xp = self.xp[3]
        self.user_xt = self.xt[3]
        self.process = process
        self.alpha, self.alpha_star = self.get_alpha(process, alpha235)
        
        # Checks are performed to ensure that only one of both streams is 
        # set to infinity.
        self.f = 0
        self.p = 0
        self.t = 0
        self.swu = 0
        self.user_f = feed
        self.user_p = product
        self.user_swu = max_swu
        
        self.n_e = n_e
        self.n_s = n_s

        self.downblend = downblend
       
        # If the desired product enrichment cannot be reached because the
        # concentration in minor isotopes is too high, then the 
        # asymptotically reached maximal enrichment is stored in this 
        # variable, else it is set to nan.
        self.maximal_enrichment = float('nan')

        return

    def check_input(self, concentrations, process, alpha235, feed, 
        product, max_swu):
        """
        Check the arguments passed when instantiating an object.
        
        This function is called upon instantiation of a Multi_isotope
        object. The arguments passed upon instantiation (or at least some
        of them) are passed to check_input, too, and they undergo some
        checks. If these are not passed corresponding errors are raised.
        """
        if alpha235<=1:
            raise ValueError("'alpha235' has to be stricly larger than 1!")
        
        try:
            concentrations['235']
        except KeyError as e:
            print('\nError:\n'
                  + 'The U235 concentrations have to be indicated as a '
                  + 'list or tuple for each stream, e.g.:\n'
                  + "'235': [0.711, 90, 0.3]")
            raise e

        for isotope in concentrations.keys():
            if isotope not in ['232', '233', '234', '235', '236']:
                msg = ('Only uranium isotopes 232, 233, 234, 235, 236, 238'
                       + 'are supported! U238 does not need to be '
                       + 'specified, its value is calculated upon '
                       + 'initialisation.\nThe keys have to be of type '
                       + "'235' and the values have to be floats (except "
                       + 'for U235 where it has to be an iterable of '
                       + 'three floats in the order "feed", "product", '
                       + '"tails".')
                raise KeyError(msg)

        if process not in ['centrifuge', 'diffusion']:
            msg = "'process' must either be 'centrifuge' or 'diffusion'!"
            raise ValueError(msg)
        
        if (max_swu==np.inf
            and feed==np.inf 
            and product==np.inf):
            msg = ("'feed', 'product' and 'max_swu' are set to infinity!\n"
                   + "At least one variable of these has to be finite!")
            raise ValueError(msg)

        if max_swu <= 0:
            raise ValueError("'max_swu' has to be stricly positive!")
        
        return
     
    def init_concentrations(self, concentrations):  
        """
        Initialise the concentration arrays.

        This function is called upon instantiation. It converts the user-
        defined concentrations into the arrays and calculates the U-238
        content in the feed. Unknown concentrations are set to 0.

        Return:
        xf -- uranium concentrations in the feed
        xp -- uranium concentrations in the product
        xt -- uranium concentrations in the tails
        """
        self.uptodate = False
        
        xf = np.zeros(6, dtype=float)
        xp = xf.copy()
        xt = xf.copy()
        
        keys = ['232', '233', '234', '235', '236']
        xf[5] = 1.

        for i, key in enumerate(keys):
            try:
                concentrations[key]
            except KeyError:
                continue

            if i != 3:
                xf[i] = concentrations[key] / 100
            else:
                xf[i] = concentrations[key][0] / 100
                xp[i] = concentrations[key][1] / 100
                xt[i] = concentrations[key][2] / 100
            xf[5] -= xf[i]
        
        if xf[5] <= 0:
            msg = ('The calculated U238 concentration in the feed is '
                   + 'negative, meaning the sum over all concentrations '
                   + 'given upon initialisation is larger than 1.')
            raise ValueError(msg)

        return xf, xp, xt
    
    def set_product_enrichment(self, xp):
        """
        Set the desired enrichment level of the U235 product in %.
        """
        if (xp < 0 
            or xp > 100
            or xp <= self.xf[3]*100):
            msg = ('"xp" has to be in the range (0, 100) and it \n'
                   + 'must be larger than the feed concentration.')
            raise ValueError(msg)
        
        self.user_xp = xp/100
        self.uptodate = False
        return True
    
    def set_tails_enrichment(self, xt):
        """
        Set the desired enrichment level of the U235 tails in %.
        """
        if (xt < 0 
            or xt > 100
            or xt >= self.xf[3]*100):
            msg = ('"xt" has to be in the range (0, 100) and it \n'
                   + 'must be smaller than the feed concentration.')
            raise ValueError(msg)
        
        self.uptodate = False
        self.user_xt = xt/100
        return True

    def set_feed_flow(self, f):
        """
        Set the desired feed flow of uranium.
        """
        if (f==np.inf
                and self.user_p==np.inf
                and self.user_swu==np.inf):
            msg = ("'f', 'user_p' and 'user_swu' are set to infinity!\n"
                   + "At least one variable of these has to be finite!")
            raise ValueError(msg)
        
        self.uptodate = False
        self.user_f = f
        return True
    
    def set_product_flow(self, p):
        """
        Set the desired product flow of uranium.
        """
        if (p==np.inf
                and self.user_f==np.inf
                and self.user_swu==np.inf):
            msg = ("'p', 'user_f' and 'user_swu' are set to infinity!\n"
                   + "At least one variable of these has to be finite!")
            raise ValueError(msg)

        self.uptodate = False
        self.user_p = p
        return True
    
    def get_alpha(self, process=None, alpha235=1.6):
        """
        Calculate the stage separation factors.
        
        Check the process used and calculate the corresponding stage 
        separation factors taking into account the cascade key weight and
        all of the isotopes. The calculations follow [1].
        """
        process = self.process if process is None else process
        if process=='centrifuge':
            # Note that the factor 1/3 has the unit of 1/atomic mass in 
            # order to keep alpha dimensionless.
            alpha = (1 + (2*self.mstar - self.isotope[3] - self.isotope) 
                         * (alpha235-1) / 3)
            alpha_star = alpha / alpha[3]**0.5
            return alpha, alpha_star
        else:
            alpha = ((2*self.mstar - self.isotope[3]) / self.isotope)**0.5
            alpha_star = alpha / alpha[3]**0.5
            return alpha, alpha_star

    def value_function(self, x):
        """
        Calculate the value of a stream with the concentrations x assuming
        that x[3] corresponds to U-235, x[5] to U-238. Based on [3].
        """
        k = (self.alpha-1) / (self.alpha[3]-1)
        abundance_ratio_235 = x[3] / x[5]
        
        if np.any(k==0.5):
            raise RuntimeError('k=0.5 in value_function() not yet implemented')
        else:
            value_func = (x / (2*k - 1)).sum()*np.log(abundance_ratio_235)
        return value_func

    def get_swu(self, f=None, p=None, t=None, xf=None, xp=None, xt=None):
        """
        Calculate the separative work units (SWU) needed in the current 
        enrichment.
        """
        f = self.f if f is None else f
        p = self.p if p is None else p 
        t = self.t if t is None else t
        xf = self.xf if xf is None else xf
        xp = self.xp if xp is None else xp
        xt = self.xt if xt is None else xt

        if ((p==np.inf or t==np.inf) and f==np.inf):
            return np.inf

        vf = self.value_function(xf)
        vp = self.value_function(xp)
        vt = self.value_function(xt)
        
        return vp*p + vt*t - vf*f
    
    def difference_concentration(self, log=False):
        """
        A helper function that returns the squared sum of relative 
        enrichment differences in product and tail.
        """
        delta_xp = (self.xp[3]-self.user_xp) / self.user_xp
        delta_xt = (self.xt[3]-self.user_xt) / self.user_xt
        
        delta = (delta_xp**2 + delta_xt**2)**0.5
        
        if log:
            msg = ("\n'difference_concentration':\n"
                   +'  delta(x_p)            {:.4e}\n'.format(delta_xp)
                   +'  delta(x_t)            {:.4e}\n'.format(delta_xt)
                   +'  delta                 {:.4e}\n'.format(delta))
            logging.debug(msg)
                
        return delta

    def calculate_concentrations(self, n_e=None, n_s=None):
        """
        Calculate the concentrations in all streams (feed, product and
        tails) for the given enrichment parameters. Following E. von Halle,
        see [2].
        """
        n_e = self.n_e if n_e is None else n_e
        n_s = self.n_s if n_s is None else n_s
        
        e = self.alpha_star**(-1) / (1 - self.alpha_star**(-n_e))
        s = self.alpha_star**(-1) / (self.alpha_star**(n_s+1) - 1)
        e_sum = (e*self.xf / (e+s)).sum()
        s_sum = (s*self.xf / (e+s)).sum()

        self.xp = e*self.xf / ((e+s) * e_sum)
        self.xt = s*self.xf / ((e+s) * s_sum)
    
        difference = self.difference_concentration()
        return difference

    def calculate_flows(self):
        """
        Calculate the material flows (feed, product and tails) for a given
        enrichment. Following E. von Halle, see [2].
        """
        e = self.alpha_star**(-1) / (1 - self.alpha_star**(-self.n_e))
        s = self.alpha_star**(-1) / (self.alpha_star**(self.n_s+1) - 1)
        e_sum = (e*self.xf / (e+s)).sum()
        s_sum = (s*self.xf / (e+s)).sum()

        p = self.user_f * e_sum
        f = self.user_p / e_sum

        if p < self.user_p:
            self.p = p
            self.f = self.user_f
        else:
            self.f = f
            self.p = self.user_p
        self.t = self.f * s_sum

        self.swu = self.get_swu()
        if self.swu > self.user_swu:
            self.swu = self.user_swu

            self.f = self.swu / (self.value_function(self.xp)*e_sum
                                 + self.value_function(self.xt)*s_sum
                                 - self.value_function(self.xf))
            self.p = self.f * e_sum
            self.t = self.f * s_sum

        return
    
    def downflow_concentration_enriching(self, stage_number):
        """
        Return the concentrations in the downflowing stream of stage 
        `stage_number` of the enriching section. This corresponds to Eq.
        (18) in [2].
        """
        self.calculate_staging()
        if (stage_number > self.n_e or stage_number < 1):
            msg = ("'stage_number' set to invalid value. It must be in the"
                   + "interval [1, {}].".format(self.n_e))
            raise ValueError(msg)
        
        bracket = 1 - self.alpha_star**(stage_number - self.n_e - 1)
        numerator = self.xp * bracket / (self.alpha_star-1)
        denominator = numerator.sum()

        return numerator / denominator
    
    def downflow_concentration_stripping(self, stage_number):
        """
        Return the concentrations in the downflowing stream of stage 
        `stage_number` of the stripping section. This corresponds to Eq.
        (34) in [2].
        """
        self.calculate_staging()
        if (stage_number > self.n_s or stage_number < 1):
            msg = ("'stage_number' set to invalid value. It must be in the"
                   + "interval [1, {}].".format(self.n_s))
            raise ValueError(msg)
        
        bracket = self.alpha_star**stage_number - 1
        numerator = self.xt * bracket / (self.alpha_star-1)
        denominator = numerator.sum()

        return numerator / denominator
    
    def interstage_enriching(self, stage_number):
        """
        Following E. von Halle Eqs. (57) and (62), see [2].
        """
        self.calculate_staging()
        
        factor1 = (self.alpha_star+1) / (self.alpha_star-1)
        factor2 = 1 - self.alpha_star**(stage_number - self.n_e - 1)
        downflowing_concentration = self.downflow_concentration_enriching(
                                                            stage_number)

        stage_flow = self.p * (self.xp*factor1*factor2).sum()
        stage_concentration_feed = (
            (self.alpha_star+1) * downflowing_concentration
            / (1 + (self.alpha_star*downflowing_concentration).sum())
        )
        
        return stage_flow, stage_concentration_feed
    
    def interstage_stripping(self, stage_number):
        """
        Following E. von Halle Eqs. (59) and (62), see [2]. However, the 
        last factor in Eq. (59) should be alpha_star^(m) - 1 instead of 
        alpha_star^(M-1).
        """
        self.calculate_staging()

        factor1 = (self.alpha_star+1) / (self.alpha_star-1)
        factor2 = self.alpha_star**stage_number - 1
        downflowing_concentration = self.downflow_concentration_stripping(
                                                            stage_number)
 
        stage_flow = self.t * (self.xt*factor1*factor2).sum()
        stage_concentration_feed = (
            (self.alpha_star+1) * stage_concentration_downflowing
            / (1 + (self.alpha_star*stage_concentration_downflowing).sum())
        )

        return stage_flow, stage_concentration_feed

    def interstage_cut(self, stage_number, section):
        """
        Return the cut of stage `stage_number`. The cut is the ratio of 
        product to feed stream. This section follows Eq. (56) in [2].
        """
        if section=="enriching":
            concentration = lambda n: self.downflow_concentration_enriching(n)
        elif section=="strippping":
            concentration = lambda n: self.downflow_concentration_stripping(n)
        else:
            raise ValueError("`section` must be either `enriching` or "
                             + "`stripping`.")
        numerator = (self.alpha_star * concentration(stage_number)).sum()
        denominator = 1 + numerator

        return numerator / denominator

    def do_downblending(self):
        EPS = 0.00005 
        USER_FEED = self.user_f
        USER_PRODUCT = self.user_p

        if (self.xp[3] - self.user_xt) < EPS:
            return
        
        blend_feed_per_product = ((self.xp[3]-self.user_xp)
                                  / (self.user_xp-self.xf[3]))
        
        if abs(self.p - self.user_p) < 1e-8:
            self.user_p /= 1 + blend_feed_per_product
            self.calculate_flows()
        elif abs(self.f - self.user_f) < 1e-8:
            e = self.alpha_star**(-1) / (1 - self.alpha_star**(-self.n_e))
            s = self.alpha_star**(-1) / (self.alpha_star**(self.n_s+1) - 1)
            e_sum = (e*self.xf / (e+s)).sum()    

            self.user_f /= 1 + blend_feed_per_product*e_sum
            self.calculate_flows()

        blend_feed = blend_feed_per_product * self.p
        if ((blend_feed + self.f > USER_FEED)
                and (abs(blend_feed + self.f - USER_FEED) > EPS)):
            e = self.alpha_star**(-1) / (1 - self.alpha_star**(-self.n_e))
            s = self.alpha_star**(-1) / (self.alpha_star**(self.n_s+1) - 1)
            e_sum = (e*self.xf / (e+s)).sum()    

            self.user_f /= 1 + blend_feed_per_product*e_sum
            self.calculate_flows()

            blend_feed = blend_feed_per_product * self.p
        
        if ((blend_feed+self.p > USER_PRODUCT)
                and (abs(blend_feed + self.p - USER_PRODUCT) > EPS)):
            self.user_p /= 1 + blend_feed_per_product
            self.calculate_flows()

            blend_feed = blend_feed_per_product * self.p

        self.xp = ((self.xp*self.p + self.xf*blend_feed) 
                   / (self.p+blend_feed))
        self.f += blend_feed
        self.p += blend_feed

        self.user_f = USER_FEED
        self.user_p = USER_PRODUCT

        return

    def calculate_staging(self, return_OptimizeResult=False):
        """
        Calculate the staging by finding the root of 
        
        sqrt((self.xp[3]-self.user_xp)**2 / self.user_xp**2
             + (self.xt[3]-self.user_xt)**2 / self.user_xt**2)
        
        using scipy's optimize.minimize function. Round the result to the 
        nearest integer value, then recalculate the cascade to get the
        final result on concentrations and streams for said integer number
        of stages. Raise RuntimeError if the optimiser does not exit
        successfully.
        """
        if self.uptodate:
            return self.n_e, self.n_s
        
        if self.process=='diffusion':
            n_init_enriching = [500, 1000, 5000]
            n_init_stripping = [100, 500, 1000, 5000]
            upper_bound = 7000
        elif self.process=='centrifuge':
            n_init_enriching = [5, 10, 50]
            n_init_stripping = [1, 5, 10, 50]
            upper_bound = 200
        else:
            msg = "'process' must either be 'centrifuge' or 'diffusion'!"
            raise ValueError(msg)
        
        concentration = lambda n=(None,None): self.calculate_concentrations(
            n_e=n[0], n_s=n[1]
        )
        lower_bound = 1
        bound = (lower_bound, upper_bound)
    
        if self.downblend:
            n = 0
            while True:
                n += 1
                self.calculate_concentrations(n, 0)
                if self.xp[3] >= self.user_xp:
                    break

            self.n_e = n
            n = 0
            while True:
                n += 1
                self.calculate_concentrations(self.n_e, n)
                if self.xt[3] <= self.user_xt:
                    break

            self.n_s = n
            self.calculate_concentrations()
            self.calculate_flows()
            self.do_downblending() 

            return

        for s in n_init_stripping:
            for e in n_init_enriching:    
                result = scipy_minimize(concentration, x0=(e, s),
                                        bounds=(bound, bound),
                                        method='L-BFGS-B',
                                        options={'gtol': 1e-15}) 
                n = result['x']
                
                msg = ("\n'calculate_staging' results:\n"
                       +"  exited successfully   {}\n".format(result['success'])
                       +"  message               {}\n".format(result['message'])
                       +"  n_e                   {:.4f}\n".format(n[0])
                       +"  n_s                   {:.4f}\n".format(n[1]))
                logging.debug(msg)
                delta = self.difference_concentration(log=True)      

                self.n_e = n[0]
                self.n_s = n[1]

                concentration()
                self.calculate_flows()
                
                if result['success'] and delta < 1e-7:
                    self.uptodate = True
                    if return_OptimizeResult:
                        return result

                    if (b'NORM_OF_PROJECTED_GRADIENT' in result['message']
                            or n[0] > 0.9*upper_bound):
                        result = scipy_minimize(concentration,
                                                x0=(10*upper_bound, s),
                                                method='L-BFGS-B',
                                                options={'gtol': 1e-15})
                        n = result['x']
                        
                        concentration(n)
                        self.calculate_flows()
                        self.maximal_enrichment = self.xp[3]
                        self.n_e = float('nan')
                        self.n_s = n[1]

                        error_msg = ('Unphysical result:\n'
                                     + 'n_enriching is larger than 0.9*'
                                     + 'upper_bound with upper_bound = '
                                     + '{}.\n'.format(upper_bound)
                                     + 'The most probable reason is that '
                                     + 'the concentration of minor isotopes '
                                     + 'is too \nhigh, making an U235 '
                                     + 'product enrichment up to the '
                                     + 'defined level impossible.\nThe '
                                     + 'maximal (asymptotical) U235 product'
                                     + 'enrichment is {:.4e} %.\n'.format(
                                            self.maximal_enrichment*100)
                                     + 'Try lowering the desired U235 '
                                     + 'enrichment below this value (e.g.,'
                                     + 'by 0.5%).')
                        raise RuntimeError(error_msg)
                    return n[0], n[1]

        error_msg = ('Optimiser did not exit successfully. Output:\n'
                     + str(result))
        raise RuntimeError(error_msg)
    
    def pprint(self, get_results=False):
        """
        Pretty print the relevant arguments passed to the object by the
        user as well as the calculated results.
        """
        self.calculate_staging() 

        print('\n--------------------------------------')
        print('Starting calculations with parameters:' )
        print('  feed               {:11.3f}'.format(self.user_f))
        print('  product            {:11.3f}'.format(self.user_p))
        print('  x_p (235)          {:11.3f}'.format(self.user_xp))
        print('  x_t (235)          {:11.3f}'.format(self.user_xt))
        print('  process            {:>11}'.format(self.process))
        print('  maximal SWU        {:11.3f}'.format(self.user_swu))
        if self.process=='centrifuge':
            print('  alpha235           {:11.3f}'.format(self.alpha[3]))
        print('\nUsed:')
        print('  feed               {:11.3f}'.format(self.f))
        print('  SWU                {:11.3f}'.format(self.swu))
        print('  enriching stages   {:11.2f}'.format(self.n_e))
        print('  stripping stages   {:11.2f}'.format(self.n_s))
        print('\nProduced:')
        print('  product            {:11.3f}'.format(self.p))
        print('  tails              {:11.3f}'.format(self.t)) 
        print('\nCompositions [%]:')
        print('  U-isotope        232         233         234'
              + '         235         236         238')
        print('  x_f     ' + (6*'{:12.4e}').format(
              self.xf[0] * 100, self.xf[1] * 100, self.xf[2] * 100, 
              self.xf[3] * 100, self.xf[4] * 100, self.xf[5] * 100))
        print('  x_p     ' + (6*'{:12.4e}').format(
              self.xp[0] * 100, self.xp[1] * 100, self.xp[2] * 100, 
              self.xp[3] * 100, self.xp[4] * 100, self.xp[5] * 100))
        print('  x_t     ' + (6*'{:12.4e}').format(
              self.xt[0] * 100, self.xt[1] * 100, self.xt[2] * 100, 
              self.xt[3] * 100, self.xt[4] * 100, self.xt[5] * 100))
        print('\n  alpha   ' + (6*'{:12.6f}').format(
              self.alpha[0], self.alpha[1], self.alpha[2], 
              self.alpha[3], self.alpha[4], self.alpha[5]))
  
        if get_results:
            return self.get_results()
        
        return

