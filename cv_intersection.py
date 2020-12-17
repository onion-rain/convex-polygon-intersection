import os
import json
import cv2
import numpy as np


def DrawPolygon(ImgShape, Polygon, Color):
    Img = np.zeros(ImgShape, np.uint8)
    try:
        cv2.fillPoly(Img, np.array(Polygon), Color)
    except:
        try:
            cv2.fillConvexPoly(Img, np.array(Polygon), Color)
        except:
            print('cant fill\n')
    return Img
 
 
def Get2PolygonIntersectArea(ImgShape, Polygon1, Polygon2):
    Img1 = DrawPolygon(ImgShape[:-1], Polygon1, 122)  # 多边形1区域填充为122
    Img2 = DrawPolygon(ImgShape[:-1], Polygon2, 133)  # 多边形2区域填充为133
    Img = Img1 + Img2
    _, OverlapImg = cv2.threshold(Img, 200, 255, cv2.THRESH_BINARY) # 根据上面的填充值，因此新图像中的像素值为255就为重叠地方
    IntersectArea = np.sum(np.greater(OverlapImg, 0)) # 求取两个多边形交叠区域面积
    return IntersectArea, OverlapImg
 
 
if __name__ == '__main__':

    # polygon1 = [[1.25,0], [0.1,1], [0.25, 1.25], [1,1]]
    # polygon2 = [[-1,1.25], [0.5,1], [-1,0], [1, 0.5], [0,2]]
    offset = 100
    scale = 50
    polygon1 = [[1*scale+offset,1*scale+offset], [-1*scale+offset,1*scale+offset], [-1*scale+offset, -1*scale+offset], [1*scale+offset,-1*scale+offset]]
    polygon2 = [[0*scale+offset,-1*scale+offset], [0*scale+offset,2*scale+offset], [2*scale+offset, 2*scale+offset], [2*scale+offset, 0*scale+offset]]

    x_coords = [p[0] for p in polygon1+polygon2]
    y_coords = [p[1] for p in polygon1+polygon2]
    background_H = (max(x_coords) - min(x_coords) + 1) + 2*offset
    background_W = (max(y_coords) - min(y_coords) + 1) + 2*offset

    background_Shape = (background_H, background_W, 3)
    IntersectArea, OverlapImg = Get2PolygonIntersectArea(background_Shape, polygon1, polygon2)
    print('IntersectArea = {}\n'.format(IntersectArea))

    cv2.imshow('OverlapImg', OverlapImg)
    Img1 = DrawPolygon(background_Shape, polygon1, (255, 0, 0))
    Img2 = DrawPolygon(background_Shape, polygon2, (0, 255, 0))
    cv2.imshow('ColorPolygons', Img1 + Img2)
    
    # # 下面使用opencv自带的函数求取一下，作为对比
    # contours, _ = cv2.findContours(OverlapImg, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    # contourArea = cv2.contourArea(contours[0])
    # print('contourArea = {}\n'.format(contourArea)) # 内部面积
    # perimeter = cv2.arcLength(contours[0], True)
    # print('contourPerimeter = {}\n'.format(perimeter)) # 轮廓面积(1像素宽的线)
    # RealContourArea = contourArea + perimeter
    # print('RealContourArea = {}\n'.format(RealContourArea))

    cv2.waitKey(0)