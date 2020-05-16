import numpy as np
import cv2
from matplotlib import pyplot as plt
import os
import math


for img_id in range(1,99):
    path="dataset/origin_99/"+str(img_id)+".jpg"
    src_img = cv2.imread(path)

    rowCount = []
    colCount = []
    for x in range(src_img.shape[0]):
        rowCountDark = 0
        for y in range(src_img.shape[1]):
            if src_img[x,y,1] < 100:
                rowCountDark += 1
        rowCount.append(rowCountDark)

    for y in range(src_img.shape[1]):
        colCountDark = 0
        for x in range(src_img.shape[0]):
            if src_img[x,y,1] < 100:
                colCountDark += 1
        colCount.append(colCountDark)

    rowMean = statistics.mean(rowCount)
    rowMedian = statistics.median(rowCount)

    colMean = statistics.mean(colCount)
    colMedian = statistics.median(colCount)

    rowTH = math.ceil(rowMedian/100.)*100
    colTH = math.floor(colMedian/100.)*100

    x1=0
    x2=0
    y1=0
    y2=0

    for index in range(len(rowCount)):
        if rowCount[index] > rowTH:
            x1 = index
            break

    for index in range( len(rowCount)-1, -1, -1):
        if rowCount[index] > rowTH:
            x2 = index
            break

    for index in range(len(colCount)):
        if colCount[index] < colTH:
            y1 = index
            break

    for index in range(len(colCount)-1, -1, -1):
        if colCount[index] < colTH:
            y2 = index
            break

    crop_img = src_img[x1:x2, y1:y2]
    cv2.imwrite("dataset/original_crop/"+str(img_id)+".jpg", crop_img)
