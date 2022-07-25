import numpy as np
from scipy.fftpack import fftshift, ifftshift
from tqdm import tqdm
from PIL import Image             # Python Imaging Library
from pathlib import Path


# Faster Fourier transform functions from the pyfftw module if
# available
try:
    from pyfftw.interfaces.scipy_fftpack import fft2, ifft2
# Otherwise use the normal scipy fftpack
except ImportError:
    import warnings
    warnings.warn("""
    Module 'pyfftw' could not be imported. Try
    'pip install pyfftw' from the Terminal. Default used,
    'fftpack' module for 2D Fourier transforms.""")
    from scipy.fftpack import fft2, ifft2

# Balanced Fast Fourier Transform.
def FFT2(I):
    [Ny,Nx]=np.shape(I)
    FT=fftshift(fft2(ifftshift(I)))/np.sqrt(Ny)/np.sqrt(Nx)
    return FT

# Balanced Inverse Fast Fourier Transform.
def iFFT2(FT):
    [Ny,Nx]=np.shape(FT)
    I=fftshift(ifft2(ifftshift(FT)))*np.sqrt(Ny)*np.sqrt(Nx)
    return I

#  Converts the matrix I to the intensity array F with values in the range [0,1].
def mat2level(I):
    I=np.array(I)
    I_Max = np.amax(I)
    I_Min = np.amin(I)
    Dim=np.shape(I)
    Ones = np.ones(Dim,dtype=np.float64)
    F = (I-Ones*I_Min)/(I_Max-I_Min)
    return(F)

# Generate SF
def speckle_gen(fsp,Ny,Nx):
    '''Function to generate an speckle field
    input = {fsp, Ny, Nx}. 
        fsp: is the pixels-radius of the frequency filter.
    output = speckle field U
    Date: 06/10/2019
    Authors: A. Federico - M. Yommi 
    '''

    # Arbitrary randon phases  
    np.random.seed(22) 
    fi = np.random.uniform(-np.pi,np.pi,(Ny,Nx))
    U = np.power(np.e,1j*fi)
    IF0 = FFT2(U)
    ny=np.array(np.linspace(-Ny/2,Ny/2,Ny))
    nx=np.array(np.linspace(-Nx/2,Nx/2,Nx))
    [x,y] = np.meshgrid(nx,ny)
    a=np.sqrt(np.multiply(y,y)+np.multiply(x,x))
    b=np.where(a>fsp,0,a)              # lowpass filter
    U = iFFT2(np.multiply(IF0,b))      # generated speckle field
    return(U)

# Iterative algorithm for the Phase Retrieval procedure (Algorithm 1, see Letter for description)
def PR(Ei,Et,thresh,b,t,NLx,NLy,steps):
    ''' Iterative phase retrieval algorithm based on alternating projections. Algorithm 1.
    Find tau for input arguments: Ei,Et,thresh,gamma,t,NLx,NLy
    input = Ei,Et,thresh,gamma,t,NLx,NLy
    outputs = tau (phase mask), RMS values: rms_angle (phase),rms_fourier(intensity)
    Date: 10/30/2019
    Authors: A. Federico - M. Yommi'''

    # Progress bar implementation
    loop=tqdm(total=steps,position=0)

    T = iFFT2(Et)
    Dim=np.shape(T)
    
    rms_fourier=[]
    rms_angle=[]
    np.random.seed(33)
    phi = np.random.uniform(-np.pi,np.pi,Dim) 
    Ai = np.abs(Ei)

    for k in range(0,steps,1):
        loop.set_description("Processing ...",format(k))
        loop.update(1)

        c = T*np.power(np.e,-1j*phi)+np.conj(T)*np.power(np.e,1j*phi)                      
        B0 = (-c+np.sqrt(np.multiply(c,c)-4*(np.multiply(np.abs(T),np.abs(T))-np.multiply(Ai,Ai))))/2 # positive solution
        B = iFFT2(FFT2(T+B0*np.power(np.e,1j*phi))*b)
        phi=np.arctan2(np.imag(B),np.real(B))
        tau = (T+B0*np.power(np.e,1j*phi))*np.power(Ei,-1)
        S = FFT2(tau*Ei)

        df=np.abs(Et)-np.abs(S)
        df=np.multiply(df,df)*t
        rms_f=np.log10(df.sum()/NLy/NLx)

        da=np.arctan2(np.imag(S),np.real(S))-np.arctan2(np.imag(Et),np.real(Et))
        da=np.multiply(da,da)*t
        rms_a=np.log10(da.sum()/NLy/NLx)
        
        rms_fourier.append(rms_f)         
        rms_angle.append(rms_a)

        if np.abs(rms_a) > thresh:
            steps=k+1
            break  

    loop.close()
    return(tau,rms_angle,rms_fourier,steps)

# Define T function
def define_T():
    ''' [T,Et,circ_in,circ_outside,domain] = define_T()
    Function to define T = IFFT[Et*t] (see Letter)
    input = none
    outputs:
    Ei, incidente complex input random field
    Et, target field
    t, target domain where Et is active
    b, exterior set of t
    thresh, value to stop the phase retrieval process
    Nx,Ny, dimension of the calculated mask
    NLx,NLy, dimension of the target field
    OffsetX,OffsetY, displacement of the origin for order zero problem
    steps, maximum number of iterations
    Date: 05/30/2019
    Authors: A. Federico - M. Yommi '''
    
    # Input conditions
    # Maximum mplitude of the Incident field
    Amp_Ei = 16  # alpha value relationship

    # Size of the target field Er
    NLy = 128
    NLx = NLy

    # Size of the phase-only mask tau  
    Nx = 270
    Ny = 320

    # Target field
    OffsetX = 0  # 72; % zero order compensation
    OffsetY = 0  # 97;

    # Threshold settings
    thresh = 20 # phase retrieval comparison
        
    # Maximum iteration number
    steps = 1000
        
    # Incident random field generation (speckle field)
    np.random.seed(17)
    phi_u=np.random.uniform(-np.pi,np.pi,[Ny,Nx])
    U = speckle_gen(135,Ny,Nx)
    AS = np.amax(np.abs(U)) 
    U = np.abs(U)/AS
    Ei = Amp_Ei*U*np.power(np.e,1j*phi_u)    
        
    # Domains and field target definitions
    uno = np.ones((Ny,Nx),dtype=np.float64)
    t = np.zeros((Ny,Nx),dtype=np.float64)
    des_y = int(Ny/2-NLy/2+OffsetY)
    has_y = int(Ny/2+NLy/2+OffsetY)
    des_x = int(Nx/2-NLx/2+OffsetX)
    has_x = int(Nx/2+NLx/2+OffsetX)

    t[des_y:has_y,des_x:has_x] = np.ones((NLy,NLx),dtype=np.float64)
    b = uno-t

    # Use BeamFormingPyProject as the default folder 
    in_file_peppers = Path.cwd() / "img" / "peppers.tif"
    in_file_baboon = Path.cwd() / "img" / "baboon.tif"

    A_t = np.zeros((Ny,Nx),dtype=np.float64)
    L = Image.open(in_file_peppers)
    L1 = L.convert('F') # 'F': 32-bit floating point scale mode
    peppers = mat2level(L1)

    Phi_t = np.zeros((Ny,Nx),dtype=np.float64)
    L = Image.open(in_file_baboon)
    L2 = L.convert('F') # 'F': 32-bit floating point scale mode
    baboon = mat2level(L2)*2*np.pi-np.pi 

    # Target field Et
    YY = peppers*np.power(np.e,1j*baboon)

    A_t[des_y:has_y,des_x:has_x] = np.abs(YY)
    Phi_t[des_y:has_y,des_x:has_x] = np.arctan2(np.imag(YY),np.real(YY))
    
    Et = A_t*np.power(np.e,1j*Phi_t) # Target beam embedded in Ny by Nx pixels
        
    return(Ei,Et,t,b,thresh,Nx,Ny,NLx,NLy,OffsetX,OffsetY,steps)

