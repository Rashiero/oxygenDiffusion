from PIL import Image,ImageSequence
import skimage
import numpy as np

OK = False
trial = 0
# Raw data
data = []
fname = 'P8_Segmented'
im = Image.open("Data/"+fname+".tif")
for i in ImageSequence.Iterator(im):
    data.append(np.array(i))

del im

data = np.array(data)
data = data.astype("int8")
Nx,Ny,Nz = data.shape

# Greatest connected component
clusters = skimage.measure.label(data,connectivity=3)

freqs = np.unique(clusters,return_counts=True)
freqs = [(x,y) for x,y in zip(freqs[0],freqs[1])]
freqs.sort(key = lambda x:x[1], reverse=True)

selected_cluster = freqs[0][0] if freqs[0][0] != 0 else freqs[1][0]

clusters[clusters != selected_cluster] = 0
clusters[clusters == selected_cluster] = 1
data = clusters.copy()

# Skeleton
skeleton = skimage.morphology.skeletonize(data)
# skeleton_old = np.copy(skeleton)

while not OK:
    # Find holes
    KNN = np.zeros((Nx,Ny,Nz))
    perc_tot = Nx*Ny*Nz
    perc_count = 0
    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                for ixi in range(max(0,ix-1),min(Nx,ix+2)):
                    for iyi in range(max(0,iy-1),min(Ny,iy+2)):
                        for izi in range(max(0,iz-1),min(Nz,iz+2)):
                            if ix != ixi or iy != iyi or iz != izi:
                                KNN[ix,iy,iz]+=skeleton[ixi,iyi,izi]
                perc_count+=1
                print(f'{100*perc_count/perc_tot:.2f}%',end="\r")

    freqs = np.unique(KNN,return_counts=True)
    dfreqs = np.array(([-1.,*freqs[0]],[0,*freqs[1]]))
    dfreqs = np.array([(dfreqs[1][i+1]-dfreqs[1][i])/(dfreqs[0][i+1]-dfreqs[0][i]) for i in range(len(freqs[0]))])

    cutoff = np.where(dfreqs > 0)[0][1]
    KNN[KNN < cutoff] = 0
    KNN[KNN != 0] = 1

    KNN = skimage.morphology.dilation(KNN)

    clusters = skimage.measure.label(np.array(KNN == 0,dtype='uint8'),connectivity=3)

    freqs = np.unique(clusters,return_counts=True)
    freqs = [(x,y) for x,y in zip(freqs[0],freqs[1])]
    freqs.sort(key = lambda x:x[1], reverse=True)

    selected_cluster = freqs[0][0]
    clusters[clusters != selected_cluster] = 0

    # Filling Holes
    data = data + skimage.morphology.dilation(np.array(clusters == 0,dtype='uint8'))
    data[data!=0] = 1
#     skeleton_old = np.copy(skeleton)
    skeleton = skimage.morphology.skeletonize(data)
    trial+=1
    skimage.io.imsave(f'{fname}_{trial}.tif',skeleton)
    OK = bool(int(input("Is it OK?\n1 - True\n0 - False\n")))
#     OK = True if len(np.unique(skeleton_old-skeleton)) == 1 else 0

skimage.io.imsave(f'{fname}_skeleton.tif',skeleton)