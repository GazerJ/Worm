---
abstract: '这个程序包由两部分构成，模拟和分析，模拟部分是基于Java写的，分析部分是基于Python3写的。主要是研究里柔软的蠕虫，在形成簇的时候，簇自发旋转的一些性质.'
author:
- 姜高晓
title: Active Worm
---

UTF8gkai

背景
====

从Vicsek提出了Vicsek模型研究鸟群后，活性物质逐渐被生物物理学家关注，利用物理学方法分析细菌，生
命系统的物理学性质，如活性物质中的有效温度，超扩散，粒子数巨涨落等。

活性物质模型
------------

Vicsek
模型在速度上做了区域内平均，这一个效应可以认为是高等生命群体的群体智能。群体内部总是保持这一个运动方向，如大雁南飞。

Active Brown Particles
模型给单个的球形粒子加了一个自驱动方向，这个模型往往去模拟细菌，但是自驱动的方向是一个角度扩散，这一点往往被人诟病

Active Worm
模型（本文），一条蠕虫身上由4个球形节点构成，由弹簧链接，每个球由LJ势能作为体积排斥是，也有弯曲能。头节点存在活性自驱动速度，方向是由尾巴指向头，大小跟周围密度成反比，这个模型可以很好的描述在集体运动的时候，个体之间活性运动方向的相互影响。

程序介绍
========

程序由两部分构成，模拟部分和分析部分,下面分别讲解程序的组成和原理

模拟
----

模拟程序主要源码由3个，置于/src/main/java/下，array.java
为控制台主程序，np.java为文件读写的程序，worm.java是核心的计算程序。

array.java:每一个worm程序类指的是一次完整独立的实验，在研究问题时需要改变各种参数，进行多次实验，这里直接在array控制即可

np.java:文件都读写。

worm.java:核心类：worm，初始化参数为扩散系数，弯曲系数，活性自驱动速度。

分析
----

这是模拟

程序使用
========

Subsection title
----------------

这是总结

12345 Vicsek T Vicsek T Vicsek T