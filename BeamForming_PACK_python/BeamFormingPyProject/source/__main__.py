import numpy as np
from aux_tools import PR, define_T, FFT2, mat2level
import matplotlib
import matplotlib.pyplot as plt  
import matplotlib.gridspec as gridspec
import matplotlib.figure as Figure
from pathlib import Path
from SSIM_PIL import compare_ssim
from PIL import Image
import csv

# Initialization
[Ei,Et,t,b,thresh,Nx,Ny,NLx,NLy,OffsetX,OffsetY,steps] = define_T()

# Iterative algorithm for phase retrieval PR                   
[tau,rms_angle,rms_fourier,steps] = PR(Ei,Et,thresh,b,t,NLx,NLy,steps)


# Results #####################################################################################

# E is the final field projected by using the obtained phase mask tau, and the incidente field Ei.
phi_mask=np.arctan2(np.imag(tau),np.real(tau))

E=FFT2(Ei*np.power(np.e,1j*phi_mask))
des_y = int(Ny/2-NLy/2+OffsetY)
has_y = int(Ny/2+NLy/2+OffsetY)
des_x = int(Nx/2-NLx/2+OffsetX)
has_x = int(Nx/2+NLx/2+OffsetX)
E_final=E[des_y:has_y,des_x:has_x]

cm = 'jet'
fig = plt.figure(constrained_layout=True)
gs0 = gridspec.GridSpec(2, 1, figure=fig) # 1 fila, 2 columnas
gs1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs0[0,0])

ax = fig.add_subplot(gs1[0])
plt.title('Recovered Amplitude')
L=np.abs(E_final)
pcm = ax.pcolormesh(np.flipud(L), cmap=cm)    
fig.colorbar(pcm, ax=ax)
plt.axis('off')

ax = fig.add_subplot(gs1[1])
plt.title('Target Amplitude')
Lt=np.abs(Et[des_y:has_y,des_x:has_x])
pcm = ax.pcolormesh(np.flipud(Lt), cmap=cm)
fig.colorbar(pcm, ax=ax)
plt.axis('off')

ax = fig.add_subplot(gs1[2])
plt.title('Recovered Phase')
Lr=np.arctan2(np.imag(E_final),np.real(E_final))
pcm = ax.pcolormesh(np.flipud(Lr), cmap=cm)
fig.colorbar(pcm, ax=ax)
plt.axis('off')

ax = fig.add_subplot(gs1[3])
plt.title('Target Phase')
Lp=np.arctan2(np.imag(Et[des_y:has_y,des_x:has_x]),np.real(Et[des_y:has_y,des_x:has_x]))
pcm = ax.pcolormesh(np.flipud(Lp), cmap=cm)
fig.colorbar(pcm, ax=ax)
plt.axis('off')

# Convergence analysis
Iteration = np.linspace(0, steps, steps) 

gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1,0])

fig.add_subplot(gs2[0]) # In phase
plt.title('Phase Convergence')
plt.plot(Iteration, rms_angle)
plt.xlabel('Iteration Number')
plt.ylabel('RMS_phase(dB)')

fig.add_subplot(gs2[1]) # In amplitude
plt.title('Amplitude Convergence')
plt.plot(Iteration, rms_fourier)
plt.xlabel('Iteration Number')
plt.ylabel('RMS_amplitude(dB)')

plt.show()

# Comparison of two images using the structural similarity algorithm (SSIM).
# https://pypi.org/project/SSIM-PIL/

# Amplitude comparison using SSIM: 
image_target = Image.fromarray((mat2level(Lt)*255).astype(int), 'L') 
image_recovered = Image.fromarray((mat2level(np.abs(E_final))*255).astype(int), 'L')  
value = compare_ssim(image_target, image_recovered)
print('\n',"The SSIM index value for amplitude calculation is Q = ",'{:f}'.format(value))

# Phase comparison using SSIM: 
image_target = Image.fromarray((mat2level(Lp)*255).astype(int), 'L') 
image_recovered = Image.fromarray((mat2level(Lr)*255).astype(int), 'L')  
value = compare_ssim(image_target, image_recovered)
print('\n',"The SSIM index value for phase calculation is Q = ",'{:f}'.format(value),'\n')

# Save results in cvs files ##################################################################################################################
input = str(input("Save values of the obtained phase-mask and the input field Ei?  Y/[N] = "))
if input=='Y':

    out_file_phi_mask = Path.cwd() / "env" / "phi_mask.csv"
    out_file_Ei_amplitude = Path.cwd() / "env" / "Ei_amplitude.csv"
    out_file_Ei_phase = Path.cwd() / "env" / "Ei_phase.csv"

    # Save the recovered phase mask
    with open(out_file_phi_mask, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(phi_mask)

    # Save the input amplitude field Ei
    with open(out_file_Ei_amplitude, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(np.abs(Ei))

    # Save the input phase field Ei
    with open(out_file_Ei_phase, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(np.arctan2(np.imag(Ei),np.real(Ei)))
