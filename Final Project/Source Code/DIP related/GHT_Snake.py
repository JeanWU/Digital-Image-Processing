#!/usr/bin/env python
# coding: utf-8



import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict

import imageio

from skimage import feature
from skimage import measure
import cv2
from scipy.ndimage.filters import sobel
import math

from imageio import imread

import morphsnakes as ms


# Path
TAR_IMG = "dataset/original_crop/"
INPUT_IMG = TAR_IMG + "crop.jpg"
OUTPUT_CUTHIST = TAR_IMG + "crop_cut.jpg"

OUTPUT_EDGE = TAR_IMG + "edge.jpg"

OUTPUT_HOUGH = TAR_IMG + "hough.jpg"

# INPUT_BACKGROUNND = TAR_IMG + "result.jpg"

OUTPUT_BLACK = TAR_IMG + "black.jpg"
OUTPUT_MERGE = TAR_IMG + "merge.jpg"

# Good for the b/w test images used
MIN_CANNY_THRESHOLD = 10
MAX_CANNY_THRESHOLD = 50


# # Cut Histogram
for k in range(1):
    img = cv2.imread(INPUT_IMG, 0)

#     hist = np.zeros((256), dtype = "int")

    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i][j] < 70:
                img[i][j] = 0
#             else:
#                 img[i][j] = quantize(img[i][j], quant)


    cv2.imwrite(OUTPUT_CUTHIST, img)
#     print(k+1)


# # Canny
L, H = 200, 230


# 要拿原圖或切完hist的圖？
img = cv2.imread(INPUT_IMG, 0)

img = cv2.medianBlur(img, 5)
canny = cv2.Canny(img, L, H)

cv2.imwrite(TAR_IMG + "median_edge.jpg", canny)



img = cv2.threshold(canny, 127, 255, cv2.THRESH_BINARY)[1]
labels = measure.label(img, connectivity=2, background=0)

mask = np.zeros(img.shape, dtype="uint8")

for (i, label) in enumerate(np.unique(labels)):
    if label == 0:
        continue

    labelMask = np.zeros(img.shape, dtype="uint8")

    labelMask[labels == label] = 255

    numPixels = cv2.countNonZero(labelMask)

    if numPixels > 20:
        mask += labelMask

cv2.imwrite(OUTPUT_EDGE, mask)


# # Hough Transform
def gradient_orientation(image):
    '''
    Calculate the gradient orientation for edge point in the image
    '''
    dx = sobel(image, axis=0, mode='constant')
    dy = sobel(image, axis=1, mode='constant')
    gradient = np.arctan2(dy,dx) * 180 / np.pi
    
    return gradient


def build_r_table(image, origin):
    '''
    Build the R-table from the given shape image and a reference point
    '''
    edges = feature.canny(image, low_threshold=MIN_CANNY_THRESHOLD, high_threshold=MAX_CANNY_THRESHOLD)
    gradient = gradient_orientation(edges)
    
    r_table = defaultdict(list)
    for (i,j),value in np.ndenumerate(edges):
        if value:
            r_table[gradient[i,j]].append((origin[0]-i, origin[1]-j))

    return r_table


def accumulate_gradients(r_table, grayImage):
    '''
    Perform a General Hough Transform with the given image and R-table
    '''
    edges = feature.canny(grayImage, low_threshold=MIN_CANNY_THRESHOLD, 
                  high_threshold=MAX_CANNY_THRESHOLD)
    gradient = gradient_orientation(edges)
    
    accumulator = np.zeros(grayImage.shape)
    for (i,j),value in np.ndenumerate(edges):
        if value:
            for r in r_table[gradient[i,j]]:
                accum_i, accum_j = i+r[0], j+r[1]
                if accum_i < accumulator.shape[0] and accum_j < accumulator.shape[1]:
                    accumulator[int(accum_i), int(accum_j)] += 1
                    
    return accumulator


def general_hough_closure(reference_image):
    '''
    Generator function to create a closure with the reference image and origin
    at the center of the reference image
    
    Returns a function f, which takes a query image and returns the accumulator
    '''
    referencePoint = (reference_image.shape[0]/2, reference_image.shape[1]/2)
    r_table = build_r_table(reference_image, referencePoint)
    
    def f(query_image):
        return accumulate_gradients(r_table, query_image)
        
    return f



def distance(x1, y1, x2, y2):
    return math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2))




def n_max(a, n, pts):
    '''
    Return the N max elements and indices in a
    '''
    indices = a.ravel().argsort()[::-1]
    indices = (np.unravel_index(i, a.shape) for i in indices)
    
    rad = 50
#     pts = [[-1, (-1, -1)]]
    
    for i in indices:
        if n <= 0:
            break
        if pts[0] == [-1, (-1, -1)]:
#             print("In")
            pts[0] = [a[i], i]
            n -= 1
        else:
            check = 1
            for [x2, y2] in pts:
                if (distance(i[0], i[1], y2[0], y2[1]) <= rad):
                    check = 0
                    break
                
            if check:
                pts.append([a[i], i])
                n -= 1
#     print([(a[i], i) for i in indices])
    
#     return [(a[i], i) for i in indices]
    print(pts)
    return pts


# In[11]:


def test_general_hough(gh, reference_image, query, points):
    '''
    Uses a GH closure to detect shapes in an image and create nice output
    '''
    query_image = cv2.imread(query, 0)
    accumulator = gh(query_image)


    detection = np.zeros((query_image.shape[0],query_image.shape[1],3), np.uint8)
    
    for i in range(detection.shape[0]):
        for j in range(detection.shape[1]):
            detection[i][j] = [query_image[i][j], query_image[i][j], query_image[i][j]]

#     detection = cv2.imread(OUTPUT_CUTHIST)
    
    # Each transform take k points
    k = 3
    m = n_max(accumulator, k, points)
    y_points = [pt[1][0] for pt in m]
    x_points = [pt[1][1] for pt in m] 
#     plt.scatter(x_points, y_points, marker='o', color='r')

    for i in range(len(x_points)):
        cv2.circle(detection,(x_points[i],y_points[i]), 3, (0,0,255), -1)

    cv2.imwrite(OUTPUT_HOUGH, detection)
    
    return


reference_image1 = cv2.imread("Hough_ref/ref_fix_1.png", 0)
reference_image2 = cv2.imread("Hough_ref/ref_2.png", 0)
detect_s_1 = general_hough_closure(reference_image1)
detect_s_2 = general_hough_closure(reference_image2)

points = [[-1, (-1, -1)]]
test_general_hough(detect_s_1, reference_image1, OUTPUT_EDGE, points)
test_general_hough(detect_s_2, reference_image2, OUTPUT_EDGE, points)

# y_points = [pt[1][0] for pt in points]
# x_points = [pt[1][1] for pt in points] 

# for i in range(len(x_points)):
#         cv2.circle(detection,(x_points[i],y_points[i]), 3, (0,0,255), -1)


# # Snake
def visual_callback_2d(background, fig=None):
    """
    Returns a callback than can be passed as the argument `iter_callback`
    of `morphological_geodesic_active_contour` and
    `morphological_chan_vese` for visualizing the evolution
    of the levelsets. Only works for 2D images.
    
    Parameters
    ----------
    background : (M, N) array
        Image to be plotted as the background of the visual evolution.
    fig : matplotlib.figure.Figure
        Figure where results will be drawn. If not given, a new figure
        will be created.
    
    Returns
    -------
    callback : Python function
        A function that receives a levelset and updates the current plot
        accordingly. This can be passed as the `iter_callback` argument of
        `morphological_geodesic_active_contour` and
        `morphological_chan_vese`.
    
    """
#     img = cv2.imread(INPUT_IMG)
    
    # Prepare the visual environment.
    if fig is None:
        fig = plt.figure()
#     plt.figure(figsize=(611/600, 370/600), dpi=600)
    
#     fig.clf()
#     fig.set_size_inches((img.shape[0], img.shape[1]))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(background, cmap=plt.cm.gray)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)


#     ax2 = fig.add_subplot(1, 2, 2)
#     ax_u = ax2.imshow(np.zeros_like(background), vmin=0, vmax=1)
    plt.pause(0.001)

    def callback(levelset):
        
        if ax.collections:
            del ax.collections[0]
        ax.contour(levelset, [0.5], colors='r')
#         ax_u.set_data(levelset)
        fig.canvas.draw()
        plt.pause(0.001)

    return callback




def example_nodule(pt):
#     logging.info('Running: example_nodule (MorphGAC)...')
    
    # Load the image.
    # 用quan完之後的圖用會很糟
    img = imread(OUTPUT_CUTHIST) / 255.0
    
    # img2為背景圖
    img2 = imread(OUTPUT_BLACK) / 255.0
    
    # g(I)
    gimg = ms.inverse_gaussian_gradient(img, alpha=1000, sigma=5.48)
    
    # Initialization of the level-set.
    init_ls = ms.circle_level_set(img.shape, pt, 15)
    
    # Callback for visual plotting
    callback = visual_callback_2d(img2)
    
    # MorphGAC. 
    ms.morphological_geodesic_active_contour(gimg, iterations=45, 
                                             init_level_set=init_ls,
                                             smoothing=1, threshold=0.29,
                                             balloon=1, iter_callback=callback)




img = cv2.imread(INPUT_IMG)

mark = np.full(img.shape, 0, dtype = 'int')

cv2.imwrite(OUTPUT_BLACK, mark)




if os.environ.get('DISPLAY', '') == '':
    matplotlib.use('Agg')
    

for i in range(len(points)):
#     fig = plt.figure()
    example_nodule(points[i][1])
    plt.savefig(TAR_IMG + ""+str(i+1)+".jpg",bbox_inches='tight',pad_inches=0.0, dpi=600)

print("Done")


# # Merge
img = cv2.imread(TAR_IMG + "1.jpg")
r = np.zeros(img.shape, dtype = 'int')

for i in range(len(points)):
#     fig = plt.figure()
    
    img = cv2.imread(TAR_IMG + ""+str(i+1)+".jpg")
        
    
    r += img
cv2.imwrite(OUTPUT_MERGE, r)




img = cv2.imread(INPUT_IMG)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.imshow(img, cmap=plt.cm.gray)
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)

plt.savefig(TAR_IMG + "trans_size.jpg",bbox_inches='tight',pad_inches=0.0, dpi = 600)




img = cv2.imread(TAR_IMG + "trans_size.jpg")
img2 = cv2.imread(OUTPUT_MERGE)

for i in range(img.shape[0]):
    for j in range(img.shape[1]):
        if img2[i][j][2] > 128:
            img[i][j] = img2[i][j]
            
cv2.imwrite(TAR_IMG + "result.jpg", img)

