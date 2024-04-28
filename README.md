杨善赟+522072910017

Week 10：

整体的重构：我们将三个结构的代码整合到一个IsingSystems的里面，里面只需要改变它们的NN就可以计算不同方面的结果，并对三者分别写unittest以及CMake。
比如说，cubic, honeycomb, kagome 都是三维的坐标晶格（按照我们的划分）, 我们只需要改变它们的临近spin挑选即可。此外我们还可以增加它们的维度等等。unittest加入能量的验证部分。

我的贡献：重新修改honeycomb和kagome的代码和unittest部分，因为之前代码的写法和cubic不太一样，为了重构写在一起整体改了一遍，方便其他组员整合。

内容见week_10 branch里面或者小组的仓库里。
