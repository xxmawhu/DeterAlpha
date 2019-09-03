# DeterAlpha
粒子衰变的角分布的一种常见方式是
		1 + alpha * cos(theta)^2
其中的参数alpha通常通过拟合得到。通常有两种方法
* chi2拟合
* 极大似然法
chi2拟合法依赖区间的划分，同时会丧失部分信息，但是操作简单。极大似然法能最大化的
利用所有的信息。这个包通过创建一个类，并且和RooFit兼容，能够方便的使用。
类名为RooDeterAlphaPdf
## 初始化方法为
```c++
    RooDeterAlphaPdf(const char* name, 
                     const char* title, 
                     RooAbsReal& _CosTheta,
                     // 参数列表，为alpha, 使得此类易扩展
                     RooArgList& params, 
                     // 效率输入
                     const TString& PHSPDat, 
                     // MC积分事例数， 不会大于PHSPDat中的事例数
                     int num = 6e6, 
                     // 是否缓存积分结果，如果"是"，可以节约积分时间
                     int store = 1
                     ); 
```
## 数据格式
数据采用txt存储，目的是为了兼容不同的ROOT版本。
txt里存储的形式为每行一个事例，每个事例包含cos(theta)和weight。比如
```txt
0.232  1.00
0.042  0.99
-0.879 1.02
```
其中weight是事例的权重，比如修正因子等。

