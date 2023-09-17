#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  17 14:18:26 2023

@author: marcos
"""

import math
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

class CPT:
    class MAF():
        def __init__(self, ch, cycles):
            self.index = 0
            self.full = [float(0)] * ch
            self.moving = [[float(0)] * ch] * cycles
            self.cycles = cycles

        def feed(self, n):
            self.full = np.subtract(self.full, self.moving[self.index]).tolist()
            self.full = np.add(self.full, n).tolist()
            self.moving[self.index] = n

            self.index = self.index + 1
            if self.index >= self.cycles:
                self.index = self.index - self.cycles
            return self

        def get(self):
            return np.divide(self.full, self.cycles).tolist()

    class UIntegral(MAF):
        def __init__(self, ch, cycles):
            super().__init__(ch, cycles)
            self.old_sample = [float(0)] * ch
            self.val = [float(0)] * ch

        def feed(self, n):
            self.val = np.add(np.add(self.val, n), self.old_sample).tolist()
            self.old_sample = n
            super().feed(self.val)
            return self

        def get(self):
            return np.multiply( (np.subtract(self.val, super().get() ) ), 2 * np.pi / (2*self.cycles) ).tolist()

    def __init__(self, channels=3, sampling_rate=61440, freq=60):
        cycles = int(sampling_rate / freq)
        	# Stores the unbiased integral (û) of the instantaneous voltage.
        self.ui = self.UIntegral(channels, cycles)
        # Stores the moving average of the instantaneous active power.
        self.p = self.MAF(channels, cycles)
        # Stores the moving average of the instantaneous active power.
        self.w = self.MAF(channels, cycles)
        # Stores the square of the RMS ||U|| of the voltage multiplied by the samples.
        self.sq_u = self.MAF(channels, cycles)
        # Stores the square of the RMS ||Û||  of the unbiased integral of the voltage multiplied by the samples.
        self.sq_ui = self.MAF(channels, cycles)
        # Stores the square of the RMS of the active current.
        self.sq_ia = self.MAF(channels, cycles)
        # Stores the square of the RMS of the reactive current.
        self.sq_ir = self.MAF(channels, cycles)
        # Stores the square of the RMS of the void current.
        self.sq_iv = self.MAF(channels, cycles)
        # Stores the square of the RMS of the balanced active current.
        self.sq_iba = self.MAF(channels, cycles)
        # Stores the square of the RMS of the balanced active current.
        self.sq_ibr = self.MAF(channels, cycles)
        # Stores the square of the RMS of the unbalanced active current.
        self.sq_iua = self.MAF(channels, cycles)
        # Stores the square of the RMS of the unbalanced reactive current.
        self.sq_iur = self.MAF(channels, cycles)

    def validatem(self, v):
        for i in range(len(v)):
            if math.isfinite(v[i]) == False:
                v[i] = float(0)
        return v

    def validate(self, v):
        if math.isfinite(v) == False:
            v = float(0)
        return v

    def feed(self, u, i):
        # Compute the instantaeous active power per phase
        p = np.multiply(u, i).tolist()

        # Compute the average active power per phase
        Pp = self.p.feed(p).get()

        # Compute the unbiased integral of the voltage
        uui = self.ui.feed(u).get()

        # Compute the instantaeous reactive energy per phase
        w = np.multiply(uui, i).tolist()

        # Compute the average active power per phase
        W = self.w.feed(w).get()

        # Compute mean square of voltage
        sq_U = self.sq_u.feed(np.multiply(u, u) ).get()

        # Compute the RMS voltage
        U = np.sqrt(sq_U).tolist()

        # Compute the instantaeous active current per phase
        ia = np.divide(np.multiply(Pp, u), sq_U).tolist()

        # Validate ia because it is possible generation an invalid value
        ia = self.validatem(ia)

        # Compute the mean square of the unbiased integral of the voltage
        sq_Ui = self.sq_ui.feed(np.multiply(uui, uui) ).get()

        # Compute the instantaeous reactive current per phase
        ir = np.divide(np.multiply(W, uui), sq_Ui).tolist()

        # Validate ia because it is possible generation an invalid value
        ir = self.validatem(ir)

        # Compute instantaeous void current per phase.
        # FIXME - These subtractions are known to cause numeric errors when using floating point types,
        # alternative strategies are advised for these cases.
        # Void power and currents should be zero for some tests and they are not.
        iv = np.subtract(np.subtract(i, ia), ir).tolist()

        # Compute the total mean square of voltage
        sq_UT = sum(sq_U)

        # Compute the total root mean square of voltage
        UT = np.sqrt(sq_UT)

        UT = self.validate(UT)

        # Compute the mean square of the void current
        sq_iV = self.sq_iv.feed(np.multiply(iv, iv) ).get()

        # Compute overall void power.
        V = UT * np.sqrt(sum(sq_iV) )

        V = self.validate(V)

        if len(u) > 1:
            # Compute instantaeous balanced active current per phase
            iba = np.divide(np.multiply(u, sum(Pp) ), sq_UT).tolist()

            iba = self.validatem(iba)

            # Compute the total mean square of the unbiased integral of voltage
            sq_UiT = sum(sq_Ui)

            # Compute instantaeous balanced reactive current per phase
            ibr = np.divide(np.multiply(uui, sum(W) ), sq_UiT).tolist()

            ibr = self.validatem(ibr)

            # Compute instantaeous unbalanced active current per phase.
            iua = np.subtract(ia, iba).tolist()

            # Compute instantaeous unbalanced reactive current per phase.
            iur = np.subtract(ir, ibr).tolist()

            # Compute overall active power.
            P = UT * np.sqrt(sum(self.sq_iba.feed(np.multiply(iba, iba) ).get() ) )

            P = self.validate(P)

			# Compute overall reactive power.
            Q = UT * np.sqrt(sum(self.sq_ibr.feed(np.multiply(ibr, ibr) ).get() ) )

            Q = self.validate(Q)

            sq_iUa = self.sq_iua.feed(np.multiply(iua, iua) ).get()
            sq_iUr = self.sq_iur.feed(np.multiply(iur, iur) ).get()

            # Compute overall unbalance power.
            N = UT * np.sqrt(sum(np.add(sq_iUa, sq_iUr) ) )

            N = self.validate(N)

            # Compute overall apparent power.
            A = np.sqrt(P*P + Q*Q + N*N + V*V)

            A = self.validate(A)

            # Compute unbalance factor.
            uf = N / np.sqrt(P*P + Q*Q + N*N)

            uf = self.validate(uf)
        else:
            # Compute overall active power.
            P = UT * np.sqrt(sum(self.sq_ia.feed(np.multiply(ia, ia) ).get() ) )

            P = self.validate(P)

            # Compute overall reactive power.
            Q = UT * np.sqrt(sum(self.sq_ir.feed(np.multiply(ir, ir) ).get() ) )

            Q = self.validate(Q)

            # Compute overall apparent power.
            A = np.sqrt(P*P + Q*Q + V*V)

            A = self.validate(A);

        # Compute reactivity factor.
        rf = Q / np.sqrt(P*P + Q*Q)

        rf = self.validate(rf)

        # Compute non linearity factor.
        df = V / A

        df = self.validate(df)

        # Compute power factor.
        pf = P / A;

        pf = self.validate(pf);

        totals = {
            'P': P,   # Overall active power.
            'Q': Q,   # Overall reactive power.
            'N': N,   # Overall unbalance power.
            'V': V,   # Overall void power.
            'A': A,   # Overall apparent power.
            'rf': rf, # Reactivity factor.
            'uf': uf, # Unbalance factor.
            'df': df, # Non linearity factor.
            'pf': pf, # Power factor.
        }

        return \
{'p': p,     # Instantaeous active power per phase.
 'P': Pp,    # Average active power per phase.
 'w': w,     # Instantaeous reactive energy per phase.
 'W': W,     # Average reactive energy per phase.
 'u': u,     # Instantaeous voltage per phase.
 'U': U,     # RMS of the voltage.
 'i': i,     # Instantaeous current per phase.
 'ia': ia,   # Instantaeous active current per phase.
 'ir': ir,   # Instantaeous reactive current per phase.
 'iv': iv,   # Instantaeous void current per phase.
 'iba': iba, # Instantaeous balanced active current per phase.
 'ibr': ibr, # Instantaeous balanced reactive current per phase.
 'iua': iua, # Instantaeous unbalanced active current per phase.
 'iur': iur, # Instantaeous unbalanced reactive current per phase.
 'totals' : totals
}

if __name__ == '__main__':
    cpt = CPT()
    angle = 0#-math.pi/2;
    for i in range(4096):
        va = 127*math.sqrt(2)*math.sin(2*math.pi*i/1024)
        vb = 127*math.sqrt(2)*math.sin(2*math.pi*i/1024 - 2*math.pi/3)
        vc = 127*math.sqrt(2)*math.sin(2*math.pi*i/1024 + 2*math.pi/3)
        ia = 1*math.sqrt(2)*math.sin(2*math.pi*i/1024 + angle)
        ib = 1*math.sqrt(2)*math.sin(2*math.pi*i/1024 - 2*math.pi/3 + angle)
        ic = 1*math.sqrt(2)*math.sin(2*math.pi*i/1024 + 2*math.pi/3 + angle)
        res = cpt.feed([va, vb, vc], [ia, ib, ic])

    print(res)
