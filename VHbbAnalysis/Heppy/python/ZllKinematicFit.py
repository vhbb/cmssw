############################################################################################################
######## Author: Tizian Bluntchli (ETH Zurich)
######## Define Estimation Class which initializes non-event-specfic data in the __init__ call      ########
######## - the idea is to create an Estimation object at he start of a scrip and then call          ########
########  its members for each event to access certain variables, i.e. the est_mass, est_pts,       ########
########      scale_factors, etc.                                                                   ########
############################################################################################################

import ROOT
import numpy as np
from numpy.linalg import inv
import os

def cut_off_one_element(array, higgsind):
    result = [[],[]]
    no_higgs_jets = sum(higgsind)
    for i in xrange(len(array)):
        cur_array = []
        cur_higgs = []
        for j in xrange(len(array)):
            if i != j:
                cur_array.append(array[j])
                cur_higgs.append(higgsind[j])
        if no_higgs_jets == sum(cur_higgs):    
            result[0].append(cur_array)
            result[1].append(cur_higgs)
    return result

def all_subarrays(nJet, higgs_indices):
    result = [[],[]]
    if nJet > 8:
        safe_result = []
        for i in xrange(nJet):
            safe_result.append(i)

        result[0] = safe_result
        result[1] = higgs_indices
        return [result]

    iter_counter = 0
    pts_indices = np.zeros(nJet, dtype = int)
    for i in xrange(nJet):
        pts_indices[i] = i

    initial_arrays = cut_off_one_element(pts_indices, higgs_indices)
    result[0].append(initial_arrays[0])
    result[1].append(initial_arrays[1])

    while(nJet - iter_counter > (sum(higgs_indices)+1)):
        cur_higgs = []
        cur_array = []

        for i in xrange(len(result[0][iter_counter][:])):

            intermediate_result = cut_off_one_element(result[0][iter_counter][i], result[1][iter_counter][i])

            for j in xrange(len(intermediate_result[0])):
                cur_array.append(intermediate_result[0][j])
                cur_higgs.append(intermediate_result[1][j])

        result[0].append(cur_array)
        result[1].append(cur_higgs)
        iter_counter += 1

    end_result_jets = []
    end_result_higgs = []
    for i in xrange(len(result[0])):
        cur_length = len(result[0][i][:])
        for j in xrange(cur_length):
            end_result_jets.append(result[0][i][j])
            end_result_higgs.append(result[1][i][j])

    combined_end_result = []
    for i in xrange(len(end_result_jets)):
        a = [end_result_jets[i], end_result_higgs[i]]
        combined_end_result.append(a)

    unique_combined_result = []
    for i in xrange(len(combined_end_result)):    
        if combined_end_result[i] not in unique_combined_result:
            unique_combined_result.append(combined_end_result[i])
            
    return unique_combined_result

def ChiSquare(V, estimated_pt, original_pt):
    V_diag = np.diag(V).tolist()
    res = 0.0
    for idx in xrange(len(estimated_pt)):
        res += ((estimated_pt[idx] - original_pt[idx])**2)/(V_diag[idx])
    return res

def A_matrix(a,b, nJet):
    diagonal = np.zeros(nJet)
    for jets in xrange(nJet):
        diagonal[jets] = a
    A = np.diag(diagonal)
    return A

def V_matrix(regions, nJet, jet_pts, jet_etas, jet_flavours, SigmasFile):

    histo_strings = []
    for jet in xrange(len(regions)):
        string = "Sigmas_" + str(int(regions[jet])) + "_" + str(int(jet_flavours[jet]))
        histo_strings.append(string)

    jet_sigmas = np.zeros(nJet)
    for jet in xrange(len(jet_sigmas)):
        SigmasFile.cd()
        current_histo = ROOT.gDirectory.Get(histo_strings[jet])
        myfunc = current_histo.GetFunction("sigma_func")

        cur_sigma = myfunc.Eval(jet_pts[jet])
        if (cur_sigma < 0.00000001) & (cur_sigma > -0.00000001):
            cur_sigma = 0.001

        jet_sigmas[jet] = myfunc.Eval(jet_pts[jet])*jet_pts[jet]

    diagonal = np.zeros(nJet)
    for jet in xrange(nJet):
        diagonal[jet] = jet_sigmas[jet]**2

    V = np.diag(diagonal)

    return V

def L_matrix_and_R_vector(nJet, jet_pts, V_pt, V_eta, V_phi, V_mass, jet_phis, jet_mass, jet_etas):
    
    Lorentzvectors = []
    for jet in xrange(nJet):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(jet_pts[jet], jet_etas[jet], jet_phis[jet], jet_mass[jet])
        Lorentzvectors.append(v)

    lepton_vector = ROOT.TLorentzVector()
    lepton_vector.SetPtEtaPhiM(V_pt, V_eta, V_phi, V_mass)

    R = np.zeros(2)
    R[0] = -lepton_vector.Px()
    R[1] = -lepton_vector.Py()

    L = np.matrix([np.cos(jet_phis), np.sin(jet_phis)])

    return L, R

def LagrangianSolver(A, L, V,  R, jet_pts):
    A_tr = np.transpose(A)
    V_inv = inv(V)
    L_tr = np.transpose(L)

    C = np.dot(A_tr, np.dot(V_inv,A))
    C_inv = inv(C)

    LC1LT1 = np.dot(L, np.dot(C_inv, L_tr))
    LC1LT1_inv = inv(LC1LT1)
    F = C_inv - np.dot(C_inv, np.dot(L_tr, np.dot(LC1LT1_inv, np.dot(L, C_inv))))
    G = np.dot(LC1LT1_inv, np.dot(L,C_inv))
    H = -1.0*LC1LT1_inv

    return np.dot(np.dot(F,np.dot(A_tr,V_inv)), jet_pts) + np.dot(np.transpose(G),R)

class ZllKinematicFit():
    def __init__(self):

        self.Resolutions = ROOT.TFile.Open(os.environ["CMSSW_BASE"] + "/src/VHbbAnalysis/Heppy/data/tf/V21_SigmasFits_fixedNames.root", "READ")
        if self.Resolutions!=None and not self.Resolutions.IsZombie():
            print "[ZllKinematicFit]: resolutions successfully loaded"
            self.resolutions_found = True
        else:
            print "[ZllKinematicFit]: resolution file not found!"
            self.resolutions_found = False
        self.original_pts = "Not yet filled"
        self.Flavours = "Not yet filled"
        self.discr_lessthan2Higgsjets = 0
        self.discr_lessthan2jets = 0
        self.discr_noresolution = 0
        self.discr_HiggsPtEtaFlavour = 0
        
    # Returns a boolean value for whether the event should be discriminated or not, the original pts w. normalized mean & the 
    def discriminate_and_custompts(self, Jet_pt, Jet_pt_reg, Jet_eta, Jet_phi, Jet_mass, hJCidx, print_discriminating_reasons = False):

        if not self.resolutions_found:
            if print_discriminating_reasons:
                print "[ZllKinematicFit]: no resolutions"
            return True

        #There can't be less than 2 Higgs jets
        if len(hJCidx) < 2:
            if print_discriminating_reasons:
                print "[ZllKinematicFit]: less than two Higgs jets"
            self.discr_lessthan2Higgsjets += 1
            return True

        #Discard events where only one jet was produced
        if len(Jet_pt) < 2:
            if print_discriminating_reasons:
                print "[ZllKinematicFit]: less than two jets"
            self.discr_lessthan2jets += 1
            return True

        # eta0: abs(eta) in [0.0,1.0]
        # eta1: abs(eta) in [1.0,1.5]
        # eta2: abs(eta) in [1.5,2.0]
        # eta3: abs(eta) in [2.0,2.5]

        regions = np.zeros(len(Jet_pt))
        for jets in xrange(len(Jet_pt)):
            if (np.absolute(Jet_eta[jets]) > 0.0 and np.absolute(Jet_eta[jets]) < 1.0):
                regions[jets] = 0
            elif (np.absolute(Jet_eta[jets]) > 1.0 and np.absolute(Jet_eta[jets]) < 1.5):
                regions[jets] = 1
            elif (np.absolute(Jet_eta[jets]) > 1.5 and np.absolute(Jet_eta[jets]) < 2.0):
                regions[jets] = 2
            elif (np.absolute(Jet_eta[jets]) > 2.0 and np.absolute(Jet_eta[jets]) < 2.5):
                regions[jets] = 3
            else: 
                regions[jets] = 99

        if any(x == 99 for x in regions):
            if print_discriminating_reasons:
                print "[ZllKinematicFit]: jet with |eta| not in the 4 regions"
            self.discr_noresolution += 1
            return True

        Flavours = np.zeros(len(Jet_pt))
        for idx in hJCidx:
            Flavours[idx] = 5

        #We wanna use Jet_pt_reg for Higgs jets (which we assume to be B-flavoured) and regular Jet_pt for the rest. 
        custom_pts = np.zeros(len(Jet_pt))
        for i in xrange(len(Jet_pt)):
            custom_pts[i] = Jet_pt[i]
        for idx in hJCidx:
            custom_pts[idx] = Jet_pt_reg[idx]

        for jets in xrange(len(Jet_pt)):    
            histo_string = "Mean_" + str(int(regions[jets])) + "_" + str(int(Flavours[jets]))
            self.Resolutions.cd()
            current_histo = ROOT.gDirectory.Get(histo_string)
            myfunc = current_histo.GetFunction("mean_func")
            custom_pts[jets] = custom_pts[jets]/(myfunc.Eval(custom_pts[jets]))
            ROOT.SetOwnership(current_histo, False ) 
            ROOT.SetOwnership(myfunc, False ) 
            #myfunc.IsA().Destructor( myfunc )
            #current_histo.IsA().Destructor( current_histo )

        #Discard all entries w Higgs-tagged Jets that have PT < 20 or |eta| > 2.4
        higgs_bools = []
        for H_jets in xrange(len(hJCidx)):
            if (custom_pts[hJCidx[H_jets]] < 20) or (Jet_eta[hJCidx[H_jets]] > 2.4) or (Jet_eta[hJCidx[H_jets]] < -2.4):
                higgs_bools.append(True)
        if any(higgs_bools):
            if print_discriminating_reasons:
                print "[ZllKinematicFit]: Higgs jet w pt < 20 or |eta| > 2.4"
            self.discr_HiggsPtEtaFlavour += 1
            return True

        self.original_pts = custom_pts
        self.Flavours = Flavours
        self.regions = regions

        return False


    # Returns the estimated pts & an array containing the indices missing the deleted entry. If the estimation didnt return
    # an array with all positive pts, event.Jet_pt and np.arange(event.nJet) are returned
    def est_pts(self, Jet_pt, Jet_pt_reg, Jet_eta, Jet_phi, Jet_mass, hJCidx, V_pt, V_eta, V_phi, V_mass):
        
        #We wanna use Jet_pt_reg for Higgs jets (which we assume to be B-flavoured) and regular Jet_pt for the rest. 
        custom_pts = self.original_pts

        higgs_indices = []
        for i in xrange(len(Jet_pt)):
            if i in hJCidx:
                higgs_indices.append(1)
            else:
                higgs_indices.append(0)

        ##### Building matrices & solving Lagrangian system #####

        # R = [(- cos (V_phi) * V_pt),(-sin(V_phi)* V_pt)]
        # L is a 2 x nJet matrix, w. the first row containing the sines of the Jet_phis, the second row containing the cosines of the Jet_phis
        #Theta are the estimated values for the jet pts

        a = 1.0
        b = 0.0

        A = A_matrix(a,b, len(Jet_pt))
        
        if np.linalg.matrix_rank(A) != A.shape[0]:
            print "Matrix A is singular"
            self.estimated_pts = custom_pts
            self.indices = np.arange(len(Jet_pt))
            return False

        V = V_matrix(self.regions, len(Jet_pt), self.original_pts, Jet_eta, self.Flavours, self.Resolutions)
        
        if np.linalg.matrix_rank(V) !=V.shape[0]:
            print "Matrix V is singular"
            self.estimated_pts = custom_pts
            self.indices = np.arange(len(Jet_pt))
            return False

        L, R = L_matrix_and_R_vector(len(Jet_pt), self.original_pts, V_pt, V_eta, V_phi, V_mass, Jet_phi, Jet_mass, Jet_eta)

        Theta = LagrangianSolver(A, L, V, R, self.original_pts)

        #List containing all possible index combinations w nJet, nJet-1, ..., len(ev.hJCidx) jets without deleting any Higgs jets.
        all_arrays = all_subarrays(len(Jet_pt), higgs_indices)

        #Inialize first guess for best estimation
        lowest_ChiSquare = ChiSquare(V, Theta[0,:].tolist()[0], self.original_pts)
        max_prob = ROOT.TMath.Prob(lowest_ChiSquare, len(Jet_pt))
        corresponding_result = Theta
        corresponding_indices = np.zeros(len(Jet_pt), dtype = int)
        for i in xrange(len(Jet_pt)):
            corresponding_indices[i] = i

        if all_arrays != []:
            for index_arrays in all_arrays:
                iteration_regions = [self.regions[i] for i in index_arrays[0]]
                iteration_nJet = len(iteration_regions)
                iteration_pts = [custom_pts[i] for i in index_arrays[0]]
                iteration_etas = [Jet_eta[i] for i in index_arrays[0]]
                iteration_Flavours = [self.Flavours[i] for i in index_arrays[0]]
                iteration_phis = [Jet_phi[i] for i in index_arrays[0]]
                iteration_masses = [Jet_mass[i] for i in index_arrays[0]]
            
                V = V_matrix(iteration_regions, iteration_nJet, iteration_pts, iteration_etas, iteration_Flavours, self.Resolutions)
                if np.linalg.matrix_rank(V) !=V.shape[0]:
                    self.estimated_pts = custom_pts
                    self.indices = np.arange(len(Jet_pt))
                    return False

                L, R = L_matrix_and_R_vector(iteration_nJet, iteration_pts, V_pt, V_eta, V_phi, V_mass, iteration_phis, iteration_masses, iteration_etas)
                if np.linalg.matrix_rank(L) !=L.shape[0]:
                    self.estimated_pts = custom_pts
                    self.indices = np.arange(len(Jet_pt))
                    return False

                A = A_matrix(a,b, len(index_arrays[0]))
                if np.linalg.matrix_rank(A) !=A.shape[0]:
                    self.estimated_pts = custom_pts
                    self.indices = np.arange(len(Jet_pt))
                    return False

                iteration_result = LagrangianSolver(A, L, V, R, iteration_pts)
                iteration_ChiSquare = ChiSquare(V, iteration_result[0,:].tolist()[0], iteration_pts)
                iteration_prob = ROOT.TMath.Prob(iteration_ChiSquare, iteration_nJet)

                neg_value_or_no = False
                for pt in iteration_result[0,:].tolist()[0]:
                    if pt < 0:
                        neg_value_or_no = True

                if not neg_value_or_no:
                    if iteration_prob > max_prob:
                        corresponding_result = iteration_result
                        corresponding_indices = index_arrays[0]


        still_negative_value_left = False

        for pt in corresponding_result[0,:].tolist()[0]:
            if pt < 0:
                still_negative_value_left = True
        
        if still_negative_value_left:
            self.estimated_pts = custom_pts
            self.indices = np.arange(len(Jet_pt))
            return False

        else:
            self.estimated_pts = corresponding_result[0,:].tolist()[0]
            self.indices = corresponding_indices
            return True


    # Uses est_pts & custom_pts to calculate the scale factors, i.e. the first order approximation w. which the masses have to be rescaled after estimation
    def scale_factors(self):

        pts = self.estimated_pts
        indices = self.indices
        custom_pts = self.original_pts

        scale_factors = np.zeros(len(indices))

        for i in xrange(len(scale_factors)):

            scale_factors[i] = pts[i] / custom_pts[indices[i]]

        self.scalefactors = scale_factors

    # Uses est_pts & scale_factors to calculate the est_mass. If estimation failed this is simply mea_mass
    def mass(self, Jet_pt, Jet_pt_reg, Jet_eta, Jet_phi, Jet_mass, hJCidx, V_pt, V_eta, V_phi, V_mass):

        self.scale_factors()

        pts = self.estimated_pts
        indices = self.indices

        est_vector = ROOT.TLorentzVector()
        mea_vector = ROOT.TLorentzVector()
        mea_vector_reg = ROOT.TLorentzVector()
        scale_factors = self.scalefactors

        new_masses = np.zeros(len(indices))
        for i in xrange(len(new_masses)):
            new_masses[i] = Jet_mass[indices[i]]*scale_factors[i]

        higgs_indc = np.zeros(len(Jet_pt), dtype= int)
        for idx in hJCidx:
            higgs_indc[idx] = 1

        higgs_post_est = [higgs_indc[i] for i in indices]

        for i in xrange(len(indices)):
            if higgs_post_est[i] == 1:
                v = ROOT.TLorentzVector()
                v.SetPtEtaPhiM(pts[i], Jet_eta[indices[i]], Jet_phi[indices[i]], new_masses[i])
                est_vector += v

        for idx in hJCidx:
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(Jet_pt_reg[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx])
            mea_vector_reg += v

            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx])
            mea_vector += v
        
        self.measured_mass_reg = mea_vector_reg.M()
        self.measured_mass = mea_vector.M()
        self.estimated_mass = est_vector.M()


