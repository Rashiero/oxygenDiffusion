from PIL import Image,ImageSequence
import skimage
import numpy as np
import scipy.ndimage as ndi

def cleaning_process(incoming,cycle = 4):
    image = incoming.copy()
    for i in range(cycle):
        print(f'Dilation {i}',end='\r')
        image = skimage.morphology.binary_dilation(image)
    for i in range(cycle):
        print(f'Erosion {i}',end='\r')
        image = skimage.morphology.binary_erosion(image)
    image = ndi.binary_fill_holes(image)
    print('Processed Image')
    return image

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

del clusters,freqs,selected_cluster

cleaned_image = skimage.morphology.remove_small_objects(frames,min_size=100)

cleaned_image = ndi.binary_fill_holes(cleaned_image)

print("Cleaning")
cleaned_image = cleaning_process(cleaned_image)

skeleton = skimage.morphology.skeletonize(cleaned_image.astype('uint8'))

skimage.io.imsave(f'{fname}_skeleton.tif',skeleton)