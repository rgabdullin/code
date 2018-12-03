#!/usr/bin/python3
# Импортируем все необходимые библиотеки:
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import sys, struct
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm as tqdm

# Читаем данные
print('reading')
res = open('/home/ruslixag/Jupyter/ParallelComputing/cache/result.bin','rb').read()
params = res[:91]
res = res[96:]

(
    Lx,Ly,Lz,
    hx,hy,hz,
    Nx,Ny,Nz,_,
    Nt,compute_fps,T,ht,
    isPX,isPY,isPZ
) = struct.unpack('6d6i2d3?',params)

T = int(T)
u = np.array(struct.unpack('{}d'.format((Nx+1)*(Ny+1)*(Nz+1)*Nt),res))
u = u.reshape((Nt,Nx+1,Ny+1,Nz+1))

print(Nx,Ny,Nz)
print(Lx,Ly,Lz)
print(hx,hy,hz)
print(T,Nt,compute_fps,ht)
print(isPX,isPY,isPZ)
print(u.shape)

num_points = (Nx+1)*(Ny+1)*(Nz+1)

x = np.array(range(0,Nx+1))
y = np.array(range(0,Ny+1))
z = np.array(range(0,Nz+1))

X,Y,Z = np.meshgrid(x,y,z)

plt_X = X.reshape((-1,))
plt_Y = Y.reshape((-1,))
plt_Z = Z.reshape((-1,))
print('X[shape] =',plt_X.shape)

print('ok')

# Объявляем все глобальные переменные
global xrot         # Величина вращения по оси x
global yrot         # Величина вращения по оси y
global zrot
global ambient      # рассеянное освещение
global greencolor   # Цвет елочных иголок
global treecolor    # Цвет елочного стебля
global lightpos     # Положение источника освещения

frame = 0

# Процедура инициализации
def init():
    global xrot
    global yrot
    global zrot
    global cmap
    global lightpos
    global ambient

    xrot = 30.0
    yrot = -30.0
    zrot = 0.0
    cmap = plt.cm.get_cmap('bwr')
    ambient = (1.0, 1.0, 1.0, 1)        # Первые три числа цвет в формате RGB, а последнее - яркость
    lightpos = (0.0,0.0, -2.0)          # Положение источника освещения по осям xyz

    glClearColor(0.0, 0.0, 0.0, 1.0)                # Серый цвет для первоначальной закраски
    gluOrtho2D(-1.0, 1.0, -1.0, 1.0)  # Определяем границы рисования по горизонтали и вертикали
    # glRotatef(-90, 0.0, 0.0, 1.0)                   # Сместимся по оси Х на 90 градусов
    # glTranslatef(0.0,0.0,-10.0)
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient) # Определяем текущую модель освещения
    glEnable(GL_LIGHTING)                           # Включаем освещение
    glEnable(GL_LIGHT0)                             # Включаем один источник света
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos)     # Определяем положение источника света
    glEnable (GL_BLEND); 
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


# Процедура обработки специальных клавиш
def specialkeys(key, x, y):
    global xrot
    global yrot
    # Обработчики для клавиш со стрелками
    if key == GLUT_KEY_UP:      # Клавиша вверх
        xrot -= 1             # Уменьшаем угол вращения по оси Х
    if key == GLUT_KEY_DOWN:    # Клавиша вниз
        xrot += 1             # Увеличиваем угол вращения по оси Х
    if key == GLUT_KEY_LEFT:    # Клавиша влево
        yrot -= 1             # Уменьшаем угол вращения по оси Y
    if key == GLUT_KEY_RIGHT:   # Клавиша вправо
        yrot += 1             # Увеличиваем угол вращения по оси Y

    glutPostRedisplay()         # Вызываем процедуру перерисовки
    print(xrot,yrot)


# Процедура перерисовки
def draw():
    print(1); sys.stdout.flush()

    global xrot
    global yrot
    global lightpos
    global cmap
    global frame
    
    glClear(GL_COLOR_BUFFER_BIT)
    glPushMatrix()
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos)

    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [1,0,0])
    glutSolidSphere(0.1,30,30)

    glTranslatef(0.2, 0.0, 0.0)
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [0,1,0])
    glutSolidSphere(0.1,30,30)

    glTranslatef(0.0, 0.2, 0.0)
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [0,0,1])
    glutSolidSphere(0.1,30,30)

    glPopMatrix()
    glPushMatrix()

    glRotatef(xrot,1.0,0,0)
    glRotatef(yrot,0,1.0,0)
    
    glTranslatef(-0.5,0.5,-0.5)

    Cmin = -0.1
    Cmax = 0.1
    plt_alpha = 0.03
    plt_size = max(hx,hy,hz)/10

    C = u[frame,:,:,:]

    plt_C = C.reshape((-1,))
    plt_C = list(map(lambda x: cmap((x - Cmin)/(Cmax-Cmin),plt_alpha),plt_C))

    def idx(i,j,k):
        return k + (Nz+1)*(j + (Ny+1)*i)

    cur_pos = np.array([0,0,0])
    for i in tqdm(range(num_points)):
        new_pos = np.array((plt_X[i] / Nx, -plt_Y[i] / Ny, plt_Z[i] / Nz))
        pos_delta = new_pos - cur_pos
        cur_pos = new_pos
        loc_color = plt_C[i]

        glTranslatef(*pos_delta)
        glMaterialfv(GL_FRONT, GL_DIFFUSE, loc_color)
        glutSolidSphere(plt_size, 30, 30)

    glPopMatrix()                                               # Возвращаем сохраненное положение "камеры"
    glutSwapBuffers()                                           # Выводим все нарисованное в памяти на экран

    glutPostRedisplay()
    frame = (frame+1) % Nt


# Здесь начинается выполнение программы
# Использовать двойную буферизацию и цвета в формате RGB (Красный, Зеленый, Синий)
glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB)
glutInitWindowSize(1000, 800)
glutInitWindowPosition(10, 10)
glutInit(sys.argv)
glutCreateWindow(b"Window1")
glutDisplayFunc(draw)
glutSpecialFunc(specialkeys)
init()
glutMainLoop()