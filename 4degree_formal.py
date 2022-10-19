import sympy as sp
from sympy.abc import x, y, z, u, v, w
from fractions import Fraction
from sympy import *
import itertools
import pprint
import copy
from sympy import Poly
import types
import time
import sys
sys.setrecursionlimit(100000)  # 将默认的递归深度修改为3000

start = time.perf_counter()
file_path = r'D:/pythonProject2/4degree_formal.out.txt' #如果要保存计算过程，根据实际情况填写，如果无需保存过程，则可以注释掉.
doc = open(file_path, 'w')
print('\n声明关于多项式的变量')
s = input("请输入要证明的关于x, y, z, u的四次多项式，变元数可以低于次数但不能高于次数：其中乘用*表示，乘方用**表示\n")
# s = x**3*z+x**3*u+x**2*y**2-x**2*y*u-2*x**2*z**2-x**2*z*u+x**2*u**2+x*y**3-x*y**2*z-x*y**2*u-x*y*z**2+x*z**3-x*z*u**2+\
#     y**3*u+y**2*z**2-2*y**2*u**2+y*z**3-y*z**2*u-y*z*u**2+y*u**3+z**2*u**2+z*u**3 #这是本文的例3
s = expand(s)
s = simplify(s)
print(s)
p = Poly(s, x, y, z, u)
variate_num = 4

# 多项式子项系数
pprint.pprint(f"多项式各项系数：{p.coeffs()}")

# 多项式的每个子项里未知数的次数
pprint.pprint(f"每个子项里各变量的次数:{p.monoms()}")

# 声明关于多项式的变量
poly_degree = 4
# 多项式次数
pprint.pprint(f"多项式次数：{poly_degree}")
print('\n')

sub_coifs = p.coeffs()
sub_list = p.monoms()

# 每个子项的字典
esub_coifs = enumerate(sub_coifs)
esub_list = enumerate(sub_list)

# new start 主要计算过程
num = 0
num_list = []
full_coifs = {}
for sub in sub_list:
    sub = enumerate(sub)
    pl = []
    for j, key in sub:
        for i in range(0, key):
            pl.append(j+1)
    print(pl)
    print(pl, file=doc)

    results = []  # pl里的元素经过排列的结果
    for s in itertools.permutations(pl, poly_degree):
        results.append(s)
    results = list(set(results))
    print(results)
    print(results, file=doc)
    num = len(set(results))
    print(num)
    print(num, file=doc)
    num_list.append(num)

    real_coifs = []  #
    for num_i, sub_coif in zip(num_list, sub_coifs):
        # coif = sub_coif*pow(num_i, -1) #结果为小数型
        coif = Fraction(sub_coif, num_i) #结果为分数型
        real_coifs.append(coif)

        for result in results:
            full_coifs[result] = coif

    print(f"\n系数a_i={real_coifs}\n")
    print(f"\n系数a_i={real_coifs}\n", file=doc)

print(full_coifs)
print(full_coifs, file=doc)

print(f"\n*K={num_list}\n")
print(f"\n*K={num_list}\n", file=doc)

full_locations = []
for key in full_coifs:
    full_locations.append(list(key))
print(f"tensor各子项下标：{full_locations}\n")
print(f"tensor各子项下标：{full_locations}\n", file=doc)

real_full_coifs = []
for value in full_coifs.values():
    real_full_coifs.append(value)
print(f"tensor各子项下标对应系数：{real_full_coifs}\n")
print(f"tensor各子项下标对应系数：{real_full_coifs}\n", file=doc)
# print('\n')


x_1 = sp.symbols('x_1')
x_2 = sp.symbols('x_2')
x_3 = sp.symbols('x_3')
x_4 = sp.symbols('x_4')

y_1 = sp.symbols('y_1')
y_2 = sp.symbols('y_2')
y_3 = sp.symbols('y_3')
y_4 = sp.symbols('y_4')

z_1 = sp.symbols('z_1')
z_2 = sp.symbols('z_2')
z_3 = sp.symbols('z_3')
z_4 = sp.symbols('z_4')

u_1 = sp.symbols('u_1')
u_2 = sp.symbols('u_2')
u_3 = sp.symbols('u_3')
u_4 = sp.symbols('u_4')

f1 = ([x_1, y_1, z_1, u_1])
f2 = ([x_2, y_2, z_2, u_2])
f3 = ([x_3, y_3, z_3, u_3])
f4 = ([x_4, y_4, z_4, u_4])
F = [f1,
     f2,
     f3,
     f4]
print(f"F={F}")
print(f"F={F}", file=doc)


def construct_v(num_of_variate, e_degree):
    v = []
    for a in range(0, e_degree):
        v_ = []
        for b in range(0, num_of_variate):
            v_.append(0)
        v_[a] = 1
        v.append(v_)
    return v


# V = [[1, 0, 0, 0],
#      [0, 1, 0, 0],
#      [0, 0, 1, 0],
#      [0, 0, 0, 1]]

V = construct_v(4, 4)
print(f"V={V}")
print(f"V={V}", file=doc)

my_tensor = []
my_tensor_with_coif = 0


def multi(locats, xc_f):
    ans_sub = 1
    for (locat, item) in zip(locats, xc_f):
        ans_sub = ans_sub * item[locat-1]
    return ans_sub


def multi_coif(list1, list2):
    ss = 0
    for (a_i, b_i) in zip(list1, list2):
        ss = ss + a_i * b_i
    return ss


for locations in full_locations:
    my_tensor_i = multi(locations, F)
    my_tensor.append(my_tensor_i)
    my_tensor_with_coif = multi_coif(my_tensor, real_full_coifs)
print(f"my_tensor={my_tensor_with_coif}\n")
print(f"my_tensor={my_tensor_with_coif}\n", file=doc)


def compute_tensor_with_coif(vertex_i, tensor_with_coif):
    tensor_with_coif_ = copy.deepcopy(tensor_with_coif)
    for (list_1, list_2) in zip(F, vertex_i):
        for k_of_list_1, j_of_list_2 in zip(list_1, list_2):
            tensor_with_coif_ = tensor_with_coif_.subs(k_of_list_1, j_of_list_2)
    print(f"运算结果：{tensor_with_coif_}")
    print(f"运算结果：{tensor_with_coif_}", file=doc)
    print('\n')
    return tensor_with_coif_


def compute_new_simplex(vertex):
    print('进行剖分')
    print('进行剖分', file=doc)
    center = []
    all_simplexes = []
    for a, b, c, d in zip(vertex[0], vertex[1], vertex[2], vertex[3]):
        center.append(Fraction((a + b + c + d), 4))
    print(f"中心：{center}")
    print(f"中心：{center}", file=doc)

    for V_j in itertools.combinations(vertex, 3):
        fate_midpoint = []
        new_simplex = [center]
        # new_simplex.append(center)
        for a, b, c in zip(V_j[0], V_j[1], V_j[2]):
            fate_midpoint.append(Fraction((a + b + c), 3))
        new_simplex.append(fate_midpoint)
        for V_jj in itertools.combinations(V_j, 2):
            side_midpoint = []
            for a, b in zip(V_jj[0], V_jj[1]):
                side_midpoint.append(Fraction((a + b), 2))
            new_simplex.append(side_midpoint)
            for one_of_side_point in V_jj:
                new_simplex.append(one_of_side_point)
                new_simplex_ = new_simplex[:]
                all_simplexes.append(new_simplex_)
                new_simplex.pop()
            new_simplex.pop()
        print('\n')
    return all_simplexes, center


def judge_whether_continue(vertex, tensor_with_coif):
    n_simplexes, center_point = compute_new_simplex(vertex)
    global j
    global m
    print(f"一共几个小单形：{len(n_simplexes)}")
    print(f"一共几个小单形：{len(n_simplexes)}\n", file=doc)
    while len(n_simplexes) > 0:
        n_simplex = n_simplexes.pop()
        print(f"新的单形顶点：{n_simplex}")  # 要计算的单形
        print(f"新的单形顶点：{n_simplex}", file=doc)  # 要计算的单形
        for v_ii in itertools.product(n_simplex, repeat=4):
            tensor_with_coif_ = compute_tensor_with_coif(v_ii, tensor_with_coif)
            if tensor_with_coif_ < 0:
                need_divide_simplexes.append(n_simplex)
                m = j + 1
                break
            else:
                pass
    division_num.append(len(need_divide_simplexes))
    if division_num[j] == 0:
        j = m
        division_num[j] = division_num[-1]
    while len(need_divide_simplexes) > 0:
        try:
            print(f"需要再次剖分的单形一共有：{len(need_divide_simplexes)}个")
            print(f"\n需要再次剖分的单形一共有：{len(need_divide_simplexes)}个", file=doc)
            nd_simplex = need_divide_simplexes.pop(0)
            division_num[j] = division_num[j]-1
            print(f"需要再次剖分的单形：{nd_simplex}")
            print(f"需要再次剖分的单形：{nd_simplex}", file=doc)
            print(f"这是第{j+1}次剖分：")
            print(f"这是第{j + 1}次剖分：", file=doc)
            yield judge_whether_continue(nd_simplex, tensor_with_coif)
        except StopIteration:
            pass
    print(f"第{j+1}次剖分结果为正")
    print(f"第{j+1}次剖分结果为正", file=doc)


def tramp(gen, *args):
    g = gen(*args)
    while isinstance(g, types.GeneratorType):
        g = g.__next__()
    return g


# 计算并剖分
division_num = [0]
m = 1
for V_i in itertools.product(V, repeat=4):
    print(V_i)
    print(V_i, file=doc)
    my_tensor_with_coif__ = compute_tensor_with_coif(V_i, my_tensor_with_coif)
    if my_tensor_with_coif__ < 0:
        if len(division_num) < 2:
            try:
                need_divide_simplexes = []
                j = 1
                print(f"这是第{j}次剖分：")
                print(f"这是第{j}次剖分：", file=doc)
                tramp(judge_whether_continue, V, my_tensor_with_coif)
            except StopIteration:
                pass
        else:
            pass
print(f"经过{division_num}次剖分，得证该不等式正定性")
print(f"经过{division_num}次剖分，得证该不等式正定性", file=doc)

end = time.perf_counter()
print("运行耗时", end-start)
print("运行耗时", end-start, file=doc)
doc.close()
