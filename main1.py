import time
import json
import asyncio
import websockets
import numpy as np
import matplotlib.pyplot as plt
import copy
import math
import random


# G = "NULL"


class item:
    def __init__(self):
        self.time_1 = []
        # self.time_2 = []
        # self.phase_1 = int()
        # self.phase_2 = int()
        # self.remains_1 = float()
        # self.remains_2 = float()
        self.direction_0 = int()
        self.direction_1 = int()
        self.direction_2 = int()
        self.direction_3 = int()
        self.remain_0 = float()
        self.remain_1 = float()
        self.remain_2 = float()
        self.remain_3 = float()
        self.a_distance = float()
        self.b_distance = float()


tls_status = item()

'''
#clientId: 2pdHJpTcWpAhKw6A
#clientSecret: YZswwBkkBiBxh4Fam5n7nTD3CMxd5xxW
'''

timestamp = int(round(time.time() * 1000))  ##13位时间戳

topic1 = "{TrafficLight}/{SN-TL-20230410-F01}/properties/report"  # 上传属性主题1

subscribe = {"type": "sub", "topic": "/device/TrafficLight/SN-TL-20230410-F01/message/property/report", "parameter": {},
             "id": "001"}
sub = json.dumps(subscribe)


def cross(len_chrom, chrom, size_pop, bound):
    # 本函数完成交叉操作
    # p_cross                input  : 交叉概率
    # len_chrom              input  : 染色体的长度
    # chrom     input  : 染色体群
    # size_pop               input  : 种群规模
    # ret                   output : 交叉后的染色体
    k1 = 0.2
    k2 = 0.9
    k3 = 0.001
    k4 = 0.01
    eta_c = 1
    # 计算适应度
    f = []
    for j in range(size_pop):
        x = chrom[j]  # 解码
        # print("x: ", x)
        f.append(fitness_function(x))  # 染色体的适应度

    f_max = max(f)  # 适应度最大值
    f_min = min(f)  # 适应度最小值
    f_avg = np.mean(f)  # 适应度平均值
    # 每一轮for循环中，可能会进行一次交叉操作，染色体是随机选择的，交叉位置也是随机选择的，%但该轮for循环中是否进行交叉操作则由交叉概率决定（continue控制）
    for i in range(size_pop):
        # 随机选择两个染色体进行交叉
        pick_sp = np.random.random_sample((1, 2))
        while np.prod(pick_sp) == 0:
            pick_sp = np.random.random_sample((1, 2))

        index_array = (pick_sp * size_pop).tolist()[0]
        index = [math.ceil(val) for val in index_array]
        for j in range(len(index)):
            if index[j] >= len(chrom):
                index[j] = index[j] - 1

        f1 = fitness_function(chrom[index[0]])  # 个体适应度值
        f2 = fitness_function(chrom[index[1]])  # 个体适应度值
        f3 = min(f1, f2)  # 两者中小者
        if f3 <= f_avg:
            # k1 * (f_max - f3) / (f_max - f_avg)
            p_cross = ((k2 - k1) * (1 - i / size_pop)) / (
                    1 + np.exp(10 * (2 * pow((f_avg - f3) / (f_avg - f_min), (f_avg / f_min)) - 1))) + k1
        else:
            p_cross = k1

        # 交叉概率决定是否进行交叉
        pick_cr = random.random()
        while pick_cr == 0:
            pick_cr = random.random()

        if pick_cr > p_cross:
            continue

        flag = 0
        while flag == 0:
            tmp_chrom = copy.deepcopy(chrom)
            # 随机选择交叉位
            pick_ind = random.random()
            while pick_ind == 0:
                pick_ind = random.random()

            index_array = (pick_ind * len(sum(len_chrom)))
            pos = int(index_array)  # 随机选择进行交叉的位置，即选择第几个变量进行交叉，注意：两个染色体交叉的位置相同

            r = random.random()
            if r <= 0.5:
                betaq = (2 * r) ** (1.0 / (eta_c + 1.0))
            else:
                betaq = (0.5 / (1.0 - r)) ** (1.0 / (eta_c + 1.0))
                # pick_ind = random.random()  # 交叉开始
            v1 = tmp_chrom[index[0]][pos]
            v2 = tmp_chrom[index[1]][pos]
            # print("v1", v1, v2)
            ylow = min(v1, v2)
            yup = max(v1, v2)
            # tmp_chrom[index[0]][pos] = pick_ind * v2 + (1 - pick_ind) * v1
            # tmp_chrom[index[1]][pos] = pick_ind * v1 + (1 - pick_ind) * v2  # 交叉结束
            tmp_chrom[index[0]][pos] = 0.5 * ((1 + betaq) * v1 + (1 - betaq) * v2)
            tmp_chrom[index[1]][pos] = 0.5 * ((1 - betaq) * v1 + (1 + betaq) * v2)

            # print("tmp_chrom[index[0]][pos] before", tmp_chrom[index[0]][pos])
            # print("tmp_chrom[index[1]][pos] before", tmp_chrom[index[1]][pos])
            tmp_chrom[index[0]][pos] = round(tmp_chrom[index[0]][pos])
            tmp_chrom[index[1]][pos] = round(tmp_chrom[index[1]][pos])
            tmp_chrom[index[0]][pos] = min(max(tmp_chrom[index[0]][pos], ylow), yup)
            tmp_chrom[index[1]][pos] = min(max(tmp_chrom[index[1]][pos], ylow), yup)
            # print("tmp_chrom[index[0]][pos]", tmp_chrom[index[0]][pos])
            # print("tmp_chrom[index[1]][pos]", tmp_chrom[index[1]][pos])
            # print("change before", chrom[index[0]])
            # print("change before", chrom[index[1]])
            flag1 = test(tmp_chrom[index[0]])  # 检验染色体1的可行性
            flag2 = test(tmp_chrom[index[1]])  # 检验染色体2的可行性
            # print("flag1", flag1)
            # print("flag2", flag2)
            # 如果两个染色体不是都可行，则重新交叉
            if flag1 * flag2 == 0:
                flag = 0
            else:
                flag = 1
                chrom[index[0]][pos] = tmp_chrom[index[0]][pos]
                chrom[index[1]][pos] = tmp_chrom[index[1]][pos]
                # print("change later", chrom[index[0]])
                # print("change later", chrom[index[1]])

            # index_array = (pick_ind * sum(lenchrom)).tolist()[0]
            # pos = math.ceil(index_array)  # 随机选择进行交叉的位置，即选择第几个变量进行交叉，注意：两个染色体交叉的位置相同
            # pick = random.random()  # 交叉开始
            # v1 = tmp_chrom[index[0]][pos]
            # v2 = tmp_chrom[index[1]][pos]
            # tmp_chrom[index[0]][pos] = pick * v2 + (1 - pick) * v1
            # tmp_chrom[index[1]][pos] = pick * v1 + (1 - pick) * v2  # 交叉结束
            #
            # flag1 = test(lenchrom, bound, tmp_chrom[index[0]])  # 检验染色体1的可行性
            # flag2 = test(lenchrom, bound, tmp_chrom[index[1]])  # 检验染色体2的可行性
            #
            # if flag1 * flag2 == 0:
            #     flag = 1
            # else:
            #     flag = 0
            # # 如果两个染色体不是都可行，则重新交叉

    ret = copy.deepcopy(chrom)
    return ret


def Mutation(len_chrom, chrom, size_pop, num, max_gen, bound):
    chrom = copy.deepcopy(chrom)
    #   % 本函数完成变异操作
    #    % p_cross                input  : 变异概率
    #   % len_chrom              input  : 染色体长度
    #   % chrom     input  : 染色体群
    #  % size_pop               input  : 种群规模
    #   % opts                  input  : 变异方法的选择
    #    % pop                   input  : 当前种群的进化代数和最大的进化代数信息
    #    % bound                 input  : 每个个体的上届和下届
    #    % max_gen                input  ：最大迭代次数
    #    % num                   input  : 当前迭代次数
    #    % ret                   output : 变异后的染色体
    k1 = 0.2
    k2 = 0.9
    k3 = 0.01
    k4 = 0.1
    # 计算适应度
    f = []
    for j in range(size_pop):
        x = chrom[j]  # 解码
        f.append(fitness_function(x))  # 染色体的适应度

    f_max = max(f)  # 适应度最大值
    f_min = min(f)  # 适应度最小值
    f_avg = np.mean(f)  # 适应度平均值

    for i in range(size_pop):  # 每一轮for循环中，可能会进行一次变异操作，染色体是随机选择的，变异位置也是随机选择的，
        # 但该轮for循环中是否进行变异操作则由变异概率决定（continue控制）
        # 随机选择一个染色体进行变异
        pick_sp = random.random()
        while pick_sp == 0:
            pick_sp = random.random()

        index = math.ceil(pick_sp * size_pop)
        if index >= len(chrom):
            index -= 1
        f1 = fitness_function(chrom[index])  # 个体适应度值
        f3 = f1
        if f3 >= f_avg:
            pmutation = ((k4 - k3) - k3 * i / max_gen) / (
                    1 + np.exp(10 * ((f_avg / f_min) * (f_avg - f3) / (f_avg - f_min) - 1))) + k3
        else:
            pmutation = k3
        # 变异概率决定该轮循环是否进行变异
        pick_cr = random.random()
        if pick_cr > pmutation:
            continue
        flag = 0
        num = 0
        chrom1 = chrom[i]
        while flag == 0 and num <= 20:
            tmp_chrom = copy.deepcopy(chrom)
            # 变异位置
            pick_loc = random.random()
            while pick_loc == 0:
                pick_loc = random.random()

            pos = int((pick_loc * sum(len_chrom)).tolist()[0])  # 随机选择了染色体变异的位置，即选择了第pos个变量进行变异
            ylow = 10
            yup = 40
            eta_m = 1
            # y = tmp_chrom[i][pos]
            delta1 = 1.0 * (tmp_chrom[i][pos] - ylow) / (yup - ylow)
            delta2 = 1.0 * (yup - tmp_chrom[i][pos]) / (yup - ylow)
            pick = random.random()  # 变异开始
            mut_pow = 1.0 / (eta_m + 1.0)
            if pick <= 0.5:
                xy = 1.0 - delta1
                val = 2.0 * pick + (1.0 - 2.0 * pick) * (xy ** (eta_m + 1.0))
                deltaq = val ** mut_pow - 1.0
            else:
                xy = 1.0 - delta2
                val = 2.0 * (1.0 - pick) + 2.0 * (pick - 0.5) * (xy ** (eta_m + 1.0))
                deltaq = 1.0 - val ** mut_pow
            tmp_chrom[i][pos] = tmp_chrom[i][pos] + deltaq * (yup - ylow)
            tmp_chrom[i][pos] = min(yup, max(tmp_chrom[i][pos], ylow))

            # fg = pow((random.random() * (1 - num / maxgen)), 2)
            # if pick > 0.5:
            #     chrom[i][pos] = chrom[i][pos] + (bound[pos][1] - chrom[i][pos]) * fg
            # else:
            #     chrom[i][pos] = chrom[i][pos] - (chrom[i][pos] - bound[pos][1]) * fg

            flag1 = test(tmp_chrom[i])  # 检验染色体的可行性
            if flag1 == 1:
                flag = 1
                chrom[i][pos] = tmp_chrom[i][pos]

        #     num = num + 1  # 检验次数设置
        #
        # if num > 20:  # 如果大于20次，则不变异
        #     chrom[i] = chrom1

    ret = copy.deepcopy(chrom)
    return ret


def test(code):
    global wait_time
    if len(code) == 1:
        t = code[0]
    else:
        t = code
    phase_period = t[0] + t[1] + t[2] + t[3]  # 信号周期
    # time_consume = 500 m / 11 m/s
    time_consume = 45.45  # 正常到达路口需要消耗的时间
    # time_arrival = 0  # 经过计算到达路口时间

    if tls_status.direction_0 == 2:
        if tls_status.remain_0 > time_consume:
            wait_time = 0
        elif (tls_status.remain_0 <= time_consume) and (tls_status.remain_0 + t[1] > time_consume):
            wait_time = tls_status.remain_0 + t[1] - time_consume + t[2] + t[3]
        elif (tls_status.remain_0 + t[1] <= time_consume) and (tls_status.remain_0 + t[1] + t[2] > time_consume):
            wait_time = tls_status.remain_0 + t[1] + t[2] - time_consume + t[3]
        elif (tls_status.remain_0 + t[1] + t[2] <= time_consume) and (
                tls_status.remain_0 + t[1] + t[2] + t[3] > time_consume):
            wait_time = tls_status.remain_0 + t[1] + t[2] + t[3] - time_consume
        elif (tls_status.remain_0 + t[1] + t[2] + t[3] <= time_consume) and (
                tls_status.remain_0 + t[1] + t[2] + t[3] + t[0] > time_consume):
            wait_time = 0

    if tls_status.direction_1 == 2:
        if (tls_status.remain_1 <= time_consume) and (tls_status.remain_1 + t[2] > time_consume):
            wait_time = tls_status.remain_1 + t[2] - time_consume + t[3]
        elif (tls_status.remain_1 + t[2] <= time_consume) and (tls_status.remain_1 + t[2] + t[3] > time_consume):
            wait_time = tls_status.remain_0 + t[2] + t[3] - time_consume
        elif (tls_status.remain_1 + t[2] + t[3] <= time_consume) and (
                tls_status.remain_1 + t[2] + t[3] + t[0] > time_consume):
            wait_time = 0
        elif tls_status.remain_1 > time_consume:
            wait_time = tls_status.remain_1 - time_consume + t[2] + t[3]
        elif (tls_status.remain_1 + t[2] + t[3] + t[0] <= time_consume) and (
                tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] > time_consume):
            wait_time = tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] - time_consume + t[2] + t[3]
        elif (tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] <= time_consume) and (
                tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] > time_consume):
            wait_time = tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] - time_consume + t[3]
        elif (tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] <= time_consume) and (
                tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] + t[3] > time_consume):
            wait_time = tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] + t[3] - time_consume
        elif (tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] + t[3] <= time_consume) and (
                tls_status.remain_1 + t[2] + t[3] + t[0] + t[1] + t[2] + t[3] + t[0] > time_consume):
            wait_time = 0

    if tls_status.direction_2 == 2:
        if (tls_status.remain_2 <= time_consume) and (tls_status.remain_2 + t[3] > time_consume):
            wait_time = tls_status.remain_2 + t[3] - time_consume
        elif (tls_status.remain_2 + t[3] <= time_consume) and (tls_status.remain_2 + t[3] + t[0] > time_consume):
            wait_time = 0
        elif (tls_status.remain_2 + t[3] + t[0] <= time_consume) and (
                tls_status.remain_2 + t[3] + t[0] + t[1] > time_consume):
            wait_time = tls_status.remain_2 + t[3] + t[0] + t[1] - time_consume + t[2] + t[3]
        elif (tls_status.remain_2 + t[3] + t[0] + t[1] <= time_consume) and (
                tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] > time_consume):
            wait_time = tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] - time_consume + t[3]
        elif (tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] <= time_consume) and (
                tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] + t[3] > time_consume):
            wait_time = tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] + t[3] - time_consume
        elif (tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] + t[3] <= time_consume) and (
                tls_status.remain_2 + t[3] + t[0] + t[1] + t[2] + t[3] + t[0] > time_consume):
            wait_time = 0
        elif tls_status.remain_2 > time_consume:
            wait_time = tls_status.remain_2 - time_consume + t[3]

    if tls_status.direction_3 == 2:
        if tls_status.remain_3 > time_consume:
            wait_time = tls_status.remain_3 - time_consume
        elif (tls_status.remain_3 <= time_consume) and (tls_status.remain_3 + t[0] > time_consume):
            wait_time = 0
        elif (tls_status.remain_3 + t[0] <= time_consume) and (tls_status.remain_3 + t[0] + t[1] > time_consume):
            wait_time = tls_status.remain_3 + t[0] + t[1] - time_consume + t[2] + t[3]
        elif (tls_status.remain_3 + t[0] + t[1] <= time_consume) and (
                tls_status.remain_3 + t[0] + t[1] + t[2] > time_consume):
            wait_time = tls_status.remain_3 + t[0] + t[1] + t[2] - time_consume + t[3]
        elif (tls_status.remain_3 + t[0] + t[1] + t[2] <= time_consume) and (
                tls_status.remain_3 + t[0] + t[1] + t[2] + t[3] > time_consume):
            wait_time = tls_status.remain_3 + t[0] + t[1] + t[2] + t[3] - time_consume
        elif (tls_status.remain_3 + t[0] + t[1] + t[2] + t[3] <= time_consume) and (
                tls_status.remain_3 + t[0] + t[1] + t[2] + t[3] + t[0] > time_consume):
            wait_time = 0

        if wait_time < 0:
            wait_time = 100

    return wait_time

    # time_arrival = tls_status.remain_0 + t[1] + t[2] + t[3]
    # if tls_status.direction_1 == 2:
    #     time_arrival = tls_status.remain_1 + t[2] + t[3]
    # if tls_status.direction_2 == 2:
    #     time_arrival = tls_status.remain_2 + t[3]
    # if tls_status.direction_3 == 2:
    #     time_arrival = tls_status.remain_3
    # flag = 1
    # if (t[0] < bound[0][0]) or (t[0] > bound[0][1]) or (t[1] < bound[1][0]) or (t[1] > bound[1][1]) or (
    #         t[2] < bound[2][0]) or (t[2] > bound[2][1]) or (t[3] < bound[3][0]) or (t[3] > bound[3][1]) or (
    #         phase_period > 150) or (phase_period < 50) or (time_consume - time_arrival > t[0]) or (time_consume - time_arrival < 0):
    #     flag = 0
    # # if (bound[0][0] <= t[0] <= bound[0][1]) & ( bound[1][0]<= t[1] <= bound[1][1]) & ( bound[2][0]<= t[2] <= bound[2][1]) & ( bound[3][0]<= t[3] <= bound[3][1]) & (90.0 <= t[0] + t[1] +t[2] + t[3] <= 130.0) :
    # # flag = 1
    # return flag


def code(bound):
    ret = []
    pick = np.random.random_sample((1, len(bound)))
    bound = np.array(bound)
    ret = np.trunc(bound[:, 0] + (bound[:, 1] - bound[:, 0]) * pick)  # 线性插值，编码结果以实数向量存入ret中
    # flag = test(len_chrom, bound, ret)  # 检验染色体的可行性
    return ret.tolist()[0]


def select(individuals, size_pop):
    individuals = copy.deepcopy(individuals)
    # 该函数用于进行选择操作
    # 求适应度值倒数
    fitness1 = [1 / val for val in individuals["fitness"]]
    # individuals.fitness为个体适应度值

    # 个体选择概率
    sum_fitness = sum(fitness1)
    sum_f = [val / sum_fitness for val in fitness1]

    # 采用轮盘赌法选择新个体实际上就是选取适应度值小的
    index = []
    # size_pop为种群数
    for i in range(size_pop):
        pick = random.random()
        while pick == 0:
            pick = random.random()
        for j in range(size_pop):
            pick = pick - sum_f[j]
            if pick < 0:
                index.append(j)
                break

    # 新种群
    new_individuals = [individuals["chrom"][j] for j in index]
    # individuals.chrom为种群中个体
    individuals["chrom"] = new_individuals
    new_fitness = [individuals["fitness"][j] for j in index]
    individuals["fitness"] = new_fitness
    ret = copy.deepcopy(individuals)
    return ret


def fitness_function(x):
    # 城市交通信号系统参数
    # C = 80  # 信号周期
    time_phase = tls_status.time_1[0] + tls_status.time_1[1] + tls_status.time_1[2] + tls_status.time_1[3]

    sy = [0, 0, 0, 0]
    sy[0] = tls_status.time_1[0] / time_phase
    sy[1] = tls_status.time_1[1] / time_phase
    sy[2] = tls_status.time_1[2] / time_phase
    sy[3] = tls_status.time_1[3] / time_phase

    t1 = int(x[0])
    t2 = int(x[1])
    t3 = int(x[2])
    t4 = int(x[3])

    time_phase_cal = t1 + t2 + t3 + t4

    lamda = [0, 0, 0, 0]
    lamda[0] = t1 / time_phase_cal  # 为第0相位的绿信比
    lamda[1] = t2 / time_phase_cal  # 为第1相位的绿信比
    lamda[2] = t3 / time_phase_cal  # 为第2相位的绿信比
    lamda[3] = t4 / time_phase_cal  # 为第3相位的绿信比

    f = 0  # 适应度值初始化
    for i in range(4):
        f += abs(sy[i] - lamda[i])
    # 解决float division by zero 的bug
    # if f == 0:
    #     f = 10
    return f


async def cloud_communication():
    # global pub
    async with websockets.connect(
            "ws://aiot.cloudxin.cn/ailinks/messaging/b0808f7899c74750542b5b54f6db73e4") as websocket:
        await websocket.send(sub)
        # print(sub)
        greeting = await websocket.recv()
        # print(greeting)
        s1 = json.loads(greeting)
        payload = s1.get("payload")
        properties = payload.get("properties")
        # print(properties)
        switch = properties.get("Switch")
        # if switch == 1:
        print(properties)
        Phase1_P3 = properties.get("Phase1_P3")
        tls_status.time_1 = Phase1_P3.get("1_time")
        # tls_status.time_1 = Phase1_P3.get("1_time")
        # tls_status.time_2 = Phase1_P3.get("2_time")
        tls_status.direction_0 = Phase1_P3.get("direction_0")
        tls_status.direction_1 = Phase1_P3.get("direction_1")
        tls_status.direction_2 = Phase1_P3.get("direction_2")
        tls_status.direction_3 = Phase1_P3.get("direction_3")
        tls_status.remain_0 = Phase1_P3.get("remain_0")
        tls_status.remain_1 = Phase1_P3.get("remain_1")
        tls_status.remain_2 = Phase1_P3.get("remain_2")
        tls_status.remain_3 = Phase1_P3.get("remain_3")
        tls_status.a_distance = Phase1_P3.get("a_distance")
        tls_status.b_distance = Phase1_P3.get("b_distance")
        # tls_status.phase_1 = Phase1_P3.get("1_phase")
        # tls_status.phase_2 = Phase1_P3.get("2_phase")
        print(Phase1_P3)
        pub = main()
        # print("main success 111111111111111")
        # await websocket.send(pub)
        # # await websocket.send(pub)
        # print(pub)
        # print("pub success 222222222222222")
        # print(111111554444454444)


global publish


def main():
    global G1
    global G2
    max_gen = 100  # 进化代数，即迭代次数
    size_pop = 100  # 种群规模
    delta = 0.1
    bound = [0, 0, 0, 0]

    # 染色体设置
    len_chrom = np.ones((1, 4))  # t1、t2、t3
    if tls_status.direction_0 == 2:
        bound = [[15, 35], [15, 35], [15, 35], [15, 35]]
    if tls_status.direction_1 == 2:
        bound = [[15, 40], [15, 35], [15, 35], [15, 35]]  # 数据范围
    if tls_status.direction_2 == 2:
        bound = [[15, 40], [15, 35], [15, 35], [15, 35]]  # 数据范围
    if tls_status.direction_3 == 2:
        bound = [[15, 40], [15, 35], [15, 35], [15, 35]]
    # ---------------------------种群初始化------------------------------------
    individuals = {'fitness': np.zeros((1, size_pop)).tolist()[0], 'chrom': []}  # 将种群信息定义为字典
    avg_fitness = []  # 每一代种群的平均适应度
    best_fitness = []  # 每一代种群的最佳适应度
    best_chrom = []  # 适应度最好的染色体
    # 初始化种群
    for i in range(size_pop):
        # 随机产生一个种群
        individuals["chrom"].append(code(bound))
        # 编码（binary和grey的编码结果为一个实数，float的编码结果为一个实数向量）
        x = individuals["chrom"][i]
        # 计算适应度
        individuals["fitness"][i] = fitness_function(x)
        # 染色体的适应度
    # 找最好的染色体
    best_fitness = min(individuals["fitness"])
    best_index = individuals["fitness"].index(best_fitness)
    best_chrom = individuals["chrom"][best_index]
    # 最好的染色体
    # 记录每一代进化中最好的适应度和平均适应度
    trace = [best_fitness]
    # 迭代求解最佳初始阀值和权值
    # 进化开始
    for i in range(max_gen):
        print('迭代次数： {} '.format(i + 1))
        # 选择
        print("select")
        individuals = select(individuals, size_pop)
        # 交叉
        print("cross")
        individuals["chrom"] = cross(len_chrom, individuals["chrom"], size_pop, bound)
        # 变异
        print("Mutation")
        individuals["chrom"] = Mutation(len_chrom, individuals["chrom"], size_pop, i, max_gen, bound)

        # 计算适应
        for j in range(size_pop):
            x = individuals["chrom"][j]  # 解码
            individuals["fitness"][j] = fitness_function(x)  # 染色体的适应度

        f_max = max(individuals["fitness"])  # 适应度最大值
        f_min = min(individuals["fitness"])  # 适应度最小值
        f_avg = np.mean(individuals["fitness"])  # 适应度平均值
        vals = [(val + abs(f_min)) / (f_max + f_min + delta) for val in individuals["fitness"]]  # 适应度标定
        individuals["fitness"] = vals
        # 找到最小和最大适应度的染色体及它们在种群中的位置
        new_best_fitness = min(individuals["fitness"])
        new_best_index = individuals["fitness"].index(new_best_fitness)
        worst_fitness = max(individuals["fitness"])
        worst_index = individuals["fitness"].index(worst_fitness)

        # 代替上一次进化中最好的染色体
        if best_fitness > new_best_fitness:
            best_fitness = new_best_fitness
            best_chrom = individuals["chrom"][new_best_index]

        individuals["chrom"][worst_index] = best_chrom  # 剔除最差个体
        trace.append(best_fitness)  # 记录每一代进化中最好的适应度

    x = best_chrom  # 最佳个体值
    D = fitness_function(best_chrom)  # 延误误差D
    print("绿信比差D", D)
    # E = D/sum(sum(q))     # 平均延误E
    # print("平均延误E",E)
    # 遗传算法结果分析
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    # # plt.rcParams['font.sans-serif'] = ['KaiTi']   # 指定默认字体
    # plt.rcParams['axes.unicode_minus'] = False
    # fig, ax = plt.subplots(1, 1)
    # plt.plot([i for i in range(len(trace))], trace, 'b--')
    # plt.title('适应度曲线  ' '终止代数＝{}'.format(max_gen))
    # plt.xlabel('进化代数')
    # plt.ylabel('适应度')
    # plt.legend('fz最佳适应度')
    # plt.show()
    G1 = [int(x[0]), int(x[1]), int(x[2]), int(x[3])]
    print("111111111111111111111111111111111", G1)
    publish = {"type": "sub", "topic": "/device-message-sender/TrafficLight/SN-TL-20220422-F01",
               "parameter": {"messageType": "INVOKE_FUNCTION", "functionId": "testPhase2",
                             "inputs": {"2_time": G1}, "headers": {"async": False}}, "id": "002"}  # 修改属性
    pub = json.dumps(publish)
    return pub


if __name__ == '__main__':

    G1 = [0, 0, 0, 0]
    G2 = [0, 0, 0, 0]
    print("======client main begin======")
    while True:
        asyncio.get_event_loop().run_until_complete(cloud_communication())  # 开启websocket线程
        # asyncio.get_event_loop().run_forever()
        # time.sleep(3)
