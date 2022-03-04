from sympy.physics.wigner import gaunt, wigner_3j



def Dipole_Coupling(l,m,l_prime,m_prime):
    coupling = pow(-1.0, m_prime)*sqrt((2.0*l_prime+1.0)*(2.0*l_block+1.0)/2.0)*wigner_3j(l_prime,1,l,0,0,0)

    coupling *= (wigner_3j(l_prime,1,l,-1*m_prime,-1,m) - wigner_3j(l_prime,1,l,-1*m_prime,1,m))

    return coupling

if __name__=="__main__":
     