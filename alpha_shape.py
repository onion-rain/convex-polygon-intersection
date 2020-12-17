def alpha_shape_2D(data, radius):
    """
    alpha shapes 算法检测边缘
    :param x: 原始点坐标集 x轴
    :param y: 原始点坐标集 y轴
    :param radius: 圆半径
    :return: 边缘点集
    """

    x = data.x
    y = data.y
    count = len(x)
    i = 0
    temp_i = i
    edge_x = []
    edge_y = []
    edge_len = 0

    while (i < count):  # 遍历点集里的所有的点
        # 根据农机作业轨迹的特性（有一定规律的，一列一列的排列）
        # 筛选 以当前点为质心，边长为 4*radius 的正方形区域 内的点 ，计算这些点与当前点的连线是否构成边界
        i_range = np.array([t for t, v in enumerate(x < x[i] + 2 * radius) if v == True])
        i_range = i_range[[t for t, v in enumerate(data.loc[i_range, 'x'] > x[i] - 2 * radius) if v == True]]
        i_range = i_range[[t for t, v in enumerate(data.loc[i_range, 'y'] < y[i] + 2 * radius) if v == True]]
        i_range = i_range[[t for t, v in enumerate(data.loc[i_range, 'y'] > y[i] - 2 * radius) if v == True]]

        """ i_range 是可能与当前点的连线组成边界线的备选点集"""
        for k in i_range:
            # 避免重复选点（一个点不能组成三个边界线，最多只能组成两条边界）
            if edge_x.count(x[k]) > 0 and edge_y.count(y[k]) > 0:
                continue

            # 计算 当前点 与 备选点 的距离
            dis = distance((x[i], y[i]), (x[k], y[k]))

            # 因为当前点 和 备选点 要在一个圆上，所有距离不能大于圆直径, 或者 当前点 与 备选点 重合
            if dis > 2 * radius or dis == 0:
                continue

            # 当前点与备选点的连线 L 线段中点
            center_x = (x[k] + x[i]) * 0.5
            center_y = (y[k] + y[i]) * 0.5

            # L 的 方向向量
            direct_x = x[k] - x[i]
            direct_y = y[k] - y[i]

            #  L 的 法向量K 及其单位化
            nomal_x = -direct_y / ma.sqrt(pow(direct_x, 2) + pow(direct_y, 2))
            nomal_y = direct_x / ma.sqrt(pow(direct_x, 2) + pow(direct_y, 2))

            # 圆心到直线L 距离
            disOfcenter = ma.sqrt(pow(radius, 2) - 0.25 * pow(dis, 2))

            if nomal_x != 0:
                a = ma.sqrt(pow(disOfcenter, 2) / (1 + pow(nomal_y / nomal_x, 2)))
                b = a * (nomal_y / nomal_x)
            else:
                a = 0
                b = radius

            # 圆心 坐标
            cicle1_x = center_x + a
            cicle1_y = center_y + b

            cicle2_x = center_x - a
            cicle2_y = center_y - b

            # 检测圆内是否有其他点
            b1 = False
            b2 = False

            # 筛选以圆心R1为质点，边长为 2*radius 的方形区域内的点集
            inSquare = np.array([t for t, v in enumerate(x < cicle1_x + radius) if v == True])
            inSquare = inSquare[[t for t, v in enumerate(data.loc[inSquare, 'x'] > cicle1_x - radius) if v == True]]
            inSquare = inSquare[[t for t, v in enumerate(data.loc[inSquare, 'y'] < cicle1_y + radius) if v == True]]
            inSquare = inSquare[[t for t, v in enumerate(data.loc[inSquare, 'y'] > cicle1_y - radius) if v == True]]
            if len(inSquare) != 0:
                for j in inSquare:
                    if j == i or j == k:  # 点集内的点 除去当前点i 和 备选点k
                        continue
                    else:
                        d = distance((x[j], y[j]), (cicle1_x, cicle1_y))  # 计算备选点k与点集内的点的距离
                        if d <= radius:  # 圆内有点 ，跳过
                            b1 = True
                            break
            # 筛选 圆R2
            inSquare = np.array([t for t, v in enumerate(x < cicle2_x + radius) if v == True])
            inSquare = inSquare[[t for t, v in enumerate(data.loc[inSquare, 'x'] > cicle2_x - radius) if v == True]]
            inSquare = inSquare[[t for t, v in enumerate(data.loc[inSquare, 'y'] < cicle2_y + radius) if v == True]]
            inSquare = inSquare[[t for t, v in enumerate(data.loc[inSquare, 'y'] > cicle2_y - radius) if v == True]]
            if len(inSquare) != 0:
                for j in inSquare:  # 与原两个点的坐标点一样
                    if j == i or j == k or distance((x[j], y[j]), (x[i], y[i])) == 0:
                        continue
                    else:
                        d = distance((x[j], y[j]), (cicle2_x, cicle2_y))
                        if d <= radius:  # 圆内有点 ，跳过
                            b2 = True
                            break

            # 圆1 或 圆2 内没有其他的点，则备选点k是边界点
            if b1 == False or b2 == False:
                if edge_x.count(x[i]) < 1 or edge_y.count(y[i]) < 1:  # 当前点未加入边界点集
                    edge_x.append(x[i])
                    edge_y.append(y[i])
                    edge.append(i)
                if edge_x.count(x[k]) < 1 or edge_y.count(y[k]) < 1:  # 备选点k已未加入边界点集
                    edge_x.append(x[k])
                    edge_y.append(y[k])
                    edge.append(k)

                temp_k = k
                break
        if edge_len < len(edge_x):  # 跳转到新的边界点
            temp_i = i
            i = temp_k
            edge_len = len(edge_x)
        else:
            temp_i = i
            i = i + 1


    edge_x.append(edge_x[0])
    edge_y.append(edge_y[0])
   return edge_x, edge_y

def distance(a, b):
    """
    计算a,b两点的直线距离
    :param a: a点坐标（大地坐标 m）
    :param b: b点坐标（大地坐标 m）
    :return: dis 直线距离
    """
    p1 = np.array(a)
    p2 = np.array(b)

    dis = ma.sqrt(sum(np.power(p2 - p1, 2)))
    return dis