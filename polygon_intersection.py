import numpy as np
import matplotlib.pyplot as plt
import math


class Edge:
    def __init__(self, point_a, point_b):
        self._support_vector = np.array(point_a)
        self._direction_vector = np.subtract(point_b, point_a)

    def get_intersection_point(self, other):
        # 重合输出也为None
        t = self._get_intersection_parameter(other)
        return None if t is None else self._get_point(t)

    def _get_point(self, parameter):
        return self._support_vector + parameter * self._direction_vector

    def _get_intersection_parameter(self, other):
        A = np.array([-self._direction_vector, other._direction_vector]).T
        if np.linalg.matrix_rank(A) < 2: # 秩小于2，说明行列式等于0，说明不可逆
        	return None
        b = np.subtract(self._support_vector, other._support_vector)
        x = np.linalg.solve(A, b)
        return x[0] if 0 <= x[0] <= 1 and 0 <= x[1] <= 1 else None


def intersect(polygon1, polygon2, sorted=1):
    """
    
    若交集为凹多边形会有bug

    原理就是求得在对方内部的顶点和边相交得到的点，依次相连得到相交多边形

    若sorted==0，则进行逆时针对点进行排序构造凸多边形

    Example: 
        polygon1 = [[1,1], [-1,1], [-1, -1], [1,-1]]
        polygon2 = [[0,0], [0,2], [2, 2], [2, 0]]

    """
    if not sorted:
        polygon1 = _sort_vertices_anti_clockwise_and_remove_duplicates(polygon1)
        polygon2 = _sort_vertices_anti_clockwise_and_remove_duplicates(polygon2)
    polygon3 = list()
    polygon3.extend(_get_vertices_lying_in_the_other_polygon(polygon1, polygon2))
    polygon3.extend(_get_edge_intersection_points(polygon1, polygon2))
    polygon3 = _sort_vertices_anti_clockwise_and_remove_duplicates(polygon3)
    return polygon1, polygon2, polygon3


def _get_vertices_lying_in_the_other_polygon(polygon1, polygon2):
    vertices = list()
    for corner in polygon1:
        if _polygon_contains_point_2(polygon2, corner):
            vertices.append(corner)
    for corner in polygon2:
        if _polygon_contains_point_2(polygon1, corner):
            vertices.append(corner)
    return vertices


def _get_edge_intersection_points(polygon1, polygon2):
    intersection_points = list()
    for i in range(len(polygon1)):
        edge1 = Edge(polygon1[i-1], polygon1[i])
        for j in range(len(polygon2)):
            edge2 = Edge(polygon2[j-1], polygon2[j])
            intersection_point = edge1.get_intersection_point(edge2)
            if intersection_point is not None:
                intersection_points.append(intersection_point)
    return intersection_points


def _polygon_contains_point(polygon, point):
    # 夹角法：判断目标点与所有边的夹角和是否为360度，为360度则在多边形内部
    angles_sum = 0
    for i in range(len(polygon)):
        angles_sum += _get_angle_from_points(point, polygon[i], polygon[i-1])
    if np.rint(angles_sum) == 360:
        return True
    else:
        return False


def _polygon_contains_point_2(polygon, point, tolerance=1e-7):
    # 引射线法：从目标点出发引一条射线，看这条射线和多边形所有边的交点数目。如果有奇数个交点，则说明在内部，如果有偶数个交点，则说明在外部
    intersection_points = [[] for i in range(4)]
    x_coords = [p[0] for p in polygon]
    y_coords = [p[1] for p in polygon]
    edge0 = Edge([min(x_coords), point[1]], point)
    edge1 = Edge([max(x_coords), point[1]], point)
    edge2 = Edge([point[0], min(y_coords)], point)
    edge3 = Edge([point[0], max(y_coords)], point)
    for i in range(len(polygon)):  # 多边形每条边
        edge = Edge(polygon[i-1], polygon[i])
        intersection_point = []
        intersection_point.append(edge0.get_intersection_point(edge))
        intersection_point.append(edge1.get_intersection_point(edge))
        intersection_point.append(edge2.get_intersection_point(edge))
        intersection_point.append(edge3.get_intersection_point(edge))
        for i in range(len(intersection_point)):  # 四条射线与该边的交点
            if intersection_point[i] is not None:
                if len(intersection_points[i]) == 0:
                    intersection_points[i].append(intersection_point[i])
                for exist_point in intersection_points[i]:
                    diff = np.subtract(exist_point, intersection_point[i])
                    if np.linalg.norm(diff, np.inf) > tolerance:
                        intersection_points[i].append(intersection_point[i])
    for point in intersection_points:
        if len(point) % 2 == 0:
            return False
    return True


def _get_angle_from_points(point_1, point_2, point_3):
    """
    根据三点坐标计算夹角
    param point_1: 点1坐标
    param point_2: 点2坐标
    param point_3: 点3坐标
    return: 返回任意角的夹角值，这里只是返回点1的夹角
    """
    a = math.sqrt((point_2[0]-point_3[0])*(point_2[0]-point_3[0])+(point_2[1]-point_3[1])*(point_2[1]-point_3[1]))
    b = math.sqrt((point_1[0]-point_3[0])*(point_1[0]-point_3[0])+(point_1[1]-point_3[1])*(point_1[1]-point_3[1]))
    c = math.sqrt((point_1[0]-point_2[0])*(point_1[0]-point_2[0])+(point_1[1]-point_2[1])*(point_1[1]-point_2[1]))
    A = math.degrees(math.acos((a*a-b*b-c*c)/(-2*b*c)))
    # B = math.degrees(math.acos((b*b-a*a-c*c)/(-2*a*c)))
    # C = math.degrees(math.acos((c*c-a*a-b*b)/(-2*a*b)))
    return A


def _sort_vertices_anti_clockwise_and_remove_duplicates(polygon, tolerance=1e-7):
    polygon = sorted(polygon, key=lambda p: _get_angle_in_radians(_get_inner_point(polygon), p))

    def vertex_not_similar_to_previous(polygon, i):
        diff = np.subtract(polygon[i-1], polygon[i])
        return i==0 or np.linalg.norm(diff, np.inf) > tolerance

    return [p for i, p in enumerate(polygon) if vertex_not_similar_to_previous(polygon, i)]


def _get_angle_in_radians(p1, p2):
    return np.arctan2(p2[1]-p1[1], p2[0]-p1[0])


def _get_inner_point(polygon):
    # 多边形外接矩形中心点
    x_coords = [p[0] for p in polygon]
    y_coords = [p[1] for p in polygon]
    return [(np.max(x_coords)+np.min(x_coords))/2., (np.max(y_coords)+np.min(y_coords))/2.]


def plot_polygon(polygon):
    polygon = list(polygon)
    polygon.append(polygon[0])
    x,y = zip(*polygon)
    plt.plot(x,y,'o-')


if __name__ == '__main__':

    polygon1 = [[np.cos(x), np.sin(x)] for x in np.random.rand(4)*2*np.pi]
    polygon2 = [[np.cos(x), np.sin(x)] for x in np.random.rand(4)*2*np.pi]

    # polygon1 = [[1.25,0], [0.1,1], [0.25, 1.25], [1,1]]
    # polygon2 = [[-1,1.25], [0.5,1], [-1,0], [1, 0.5], [0,2]]
    # polygon1 = [[1,1], [-1,1], [-1, -1], [1,-1]]
    # polygon2 = [[0,-1], [0,2], [2, 0], [2, 2]]
    polygon1, polygon2, polygon3 = intersect(polygon1, polygon2, sorted=0)

    plot_polygon(polygon1)
    plot_polygon(polygon2)
    if len(polygon3) > 0:
        plot_polygon(polygon3)
    plt.show()
    print()
