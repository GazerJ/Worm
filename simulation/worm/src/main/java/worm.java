import java.io.IOException;
import java.util.Random;


public class worm extends Thread {
    final double L = 100;
    final double L_2=L/2;
    //final double cutr=6;//80
    final double cutr=10;//400
    final double dt = 0.005;   // 0.0001 *1000*10000 10h
    final int n = 80;
    final  double lam;
    final int nStep = 11;   //每次迭代计算nStep-1 个时间步长                             // nStep  占内存   needtime 占IO
    final String idNum;
    final double T;//极限0.01，超过这个时候会发散
    final double Ts;
    final double v0;
    final String Tstr;
    final String Tsstr;
    final String v0str;
    final int needTime =400000;// 200粒子1次1.8s//75粒子0.32s //500粒子一次24秒
    final int cuttt=100;
    long start = System.currentTimeMillis();

    public worm(int i, double T,double Ts, double lam, double v0) {
        this.idNum = String.format("%d", i);
        this.T = T * Math.sqrt(dt);    //初始化的时候直接把Ts 设为传过来的参数的Ts*sqrt(dt)
        this.Ts =Ts * Math.sqrt(dt);  // Ts * Math.sqrt(dt); 初始化的时候直接把Ts 设为传过来的参数的Ts*sqrt(dt)
        this.v0 = v0;
        this.lam=lam;

        this.Tstr = String.format("%.8f", T);
        this.Tsstr = String.format("%.8f", lam);
        this.v0str = String.format("%.1f", v0);
    }


    public void run() {
        System.out.println("id:" + this.idNum + "\tT:" + this.Tstr + "\tlam:" + this.Tsstr);

        final double[][] x1 = new double[nStep][n];
        final double[][] y1 = new double[nStep][n];
        final double[][] x2 = new double[nStep][n];
        final double[][] y2 = new double[nStep][n];
        final double[][] x3 = new double[nStep][n];
        final double[][] y3 = new double[nStep][n];
        final double[][] x4 = new double[nStep][n];
        final double[][] y4 = new double[nStep][n];

        try {
            x1[0] = np.loadtxt("./input/toheart30x1.txt", n);
            x2[0] = np.loadtxt("./input/toheart30x2.txt", n);
            x3[0] = np.loadtxt("./input/toheart30x3.txt", n);
            x4[0] = np.loadtxt("./input/toheart30x4.txt", n);
            y1[0] = np.loadtxt("./input/toheart30y1.txt", n);
            y2[0] = np.loadtxt("./input/toheart30y2.txt", n);
            y3[0] = np.loadtxt("./input/toheart30y3.txt", n);
            y4[0] = np.loadtxt("./input/toheart30y4.txt", n);
            gotest(x1, y1, x2, y2, x3, y3, x4, y4);
            System.out.println(System.currentTimeMillis() - start);
            System.out.println("end");

        } catch (IOException ex) {
            ex.printStackTrace();
        }


    }


    public void gotest(double[][] x1, double[][] y1, double[][] x2, double[][] y2, double[][] x3, double[][] y3, double[][] x4, double[][] y4) throws IOException {

        final Random rnd = new Random();

        double xiX;
        double xiY;
        //Random rnd = new Random();
        double[][] X = new double[4][n];
        double[][] Y = new double[4][n];
        double[][] dx;
        double[][] dy;
        double[][] dr;
        double[][] fljcore;
        double[][][] Tcore;
        final double[][] Fx = new double[4][n];
        final double[][] Fy = new double[4][n];
        double[][] Factive;
        for (int times = 0; times < needTime; times++) {

            for (int i = 0; i < nStep - 1; i++) {

//                X = new double[][]{x1[i], x2[i], x3[i], x4[i]};
//                Y = new double[][]{y1[i], y2[i], y3[i], y4[i]};
                X[0] = x1[i];
                X[1] = x2[i];
                X[2] = x3[i];
                X[3] = x4[i];
                Y[0] = y1[i];
                Y[1] = y2[i];
                Y[2] = y3[i];
                Y[3] = y4[i];

                dx = getDx(X);
                dy = getDy(Y);
                dr = getDr(dx, dy);
                fljcore = getFljcore(dr);
                Tcore = getTcore(dr, dx, dy, X, Y);
                Fx[0] = gexFx(fljcore, dx, 0);
                Fy[0] = gexFy(fljcore, dy, 0);
                Fx[1] = gexFx(fljcore, dx, 1);
                Fy[1] = gexFy(fljcore, dy, 1);
                Fx[2] = gexFx(fljcore, dx, 2);
                Fy[2] = gexFy(fljcore, dy, 2);
                Fx[3] = gexFx(fljcore, dx, 3);
                Fy[3] = gexFy(fljcore, dy, 3);
                Factive = getFactive(dx, dy, dr);
                for (int j = 0; j < n; j++) {
                    //xiX = T * rnd.nextGaussian(); //X方向  空间扩散
                    //xiY = T * rnd.nextGaussian(); //Y方向  空间扩散

                    x1[i + 1][j] = Math.abs((x1[i][j] + (Factive[j][0] + Tcore[0][j][0] + Fx[0][j]) * dt + 1*T * rnd.nextGaussian()+L) % L); //%是取余，也就是周期边界条件
                    y1[i + 1][j] = Math.abs((y1[i][j] + (Factive[j][1] + Tcore[0][j][1] + Fy[0][j]) * dt + 1*T * rnd.nextGaussian()+L) % L);
                    x2[i + 1][j] = Math.abs((x2[i][j] + (Tcore[1][j][0] + Fx[1][j]) * dt + 1*T * rnd.nextGaussian()+L) % L);
                    y2[i + 1][j] = Math.abs((y2[i][j] + (Tcore[1][j][1] + Fy[1][j]) * dt + 1*T*rnd.nextGaussian()+L) % L);
                    x3[i + 1][j] = Math.abs((x3[i][j] + (Tcore[2][j][0] + Fx[2][j]) * dt + 1*T * rnd.nextGaussian()+L) % L);
                    y3[i + 1][j] = Math.abs((y3[i][j] + (Tcore[2][j][1] + Fy[2][j]) * dt + 1*T * rnd.nextGaussian()+L) % L);
                    x4[i + 1][j] = Math.abs((x4[i][j] + (Tcore[3][j][0] + Fx[3][j]) * dt + 1*T * rnd.nextGaussian()+L) % L);
                    y4[i + 1][j] = Math.abs((y4[i][j] + (Tcore[3][j][1] + Fy[3][j]) * dt + 1*T * rnd.nextGaussian()+L) % L);
                }

            }
            if(times %cuttt ==0 ) {
                System.out.print("ID : ");
                System.out.print(this.idNum);
                System.out.print("times : ");
                System.out.print(times);
                System.out.print("    have run time: ");
                System.out.print(System.currentTimeMillis() - start);
                System.out.print("    left time: ");
                System.out.print((double)(System.currentTimeMillis()/1000 - start/1000)/(double)(3600*times) *(needTime-times));
                System.out.println("    hours ");


                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "x1.txt", x1, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "x2.txt", x2, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "x3.txt", x3, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "x4.txt", x4, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "y1.txt", y1, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "y2.txt", y2, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "y3.txt", y3, nStep);
                np.savetxtCutAdd("./output/data" + this.idNum + "/T" + Tstr + "Ts" + Tsstr + "v0" + v0str + "y4.txt", y4, nStep);

            }
            //System.out.print("ID\t"+this.idNum+"times\t"+String.format("%d", times)+'\t');
            //System.out.println(System.currentTimeMillis() - start);

//            double[] tmpx1 = x1[nStep - 1].clone();
//            double[] tmpx2 = x2[nStep - 1].clone();
//            double[] tmpx3 = x3[nStep - 1].clone();
//            double[] tmpx4 = x4[nStep - 1].clone();
//            double[] tmpy1 = y1[nStep - 1].clone();
//            double[] tmpy2 = y2[nStep - 1].clone();
//            double[] tmpy3 = y3[nStep - 1].clone();
//            double[] tmpy4 = y4[nStep - 1].clone();

            //np.fill(x1,x2,x3,x4,y1,y2,y3,y4);

            x1[0] = x1[nStep - 1].clone();

            x2[0] = x2[nStep - 1].clone();
            x3[0] = x3[nStep - 1].clone();
            x4[0] = x4[nStep - 1].clone();
            y1[0] = y1[nStep - 1].clone();
            y2[0] = y2[nStep - 1].clone();
            y3[0] = y3[nStep - 1].clone();
            y4[0] = y4[nStep - 1].clone();


        }
    }

    public double[][] getFactive(double[][] dx, double[][] dy, double[][] dr) {
        final double[][] fActive = new double[n][2]; //n条蠕虫活性力，2个方向，所以是n*2的数组
        final double[] p = new double[2];////  一条蠕虫的取向数组
        final int n3 = n * 3;             // 定义的常数n3就是n*3     优化程序而定义的尽量减少过程中的重复计算
        int n3i;                          // 定义的常数n3i就是n*3+i     优化程序而定义的尽量减少过程中的重复计算

        double fCore;
        double labelx = 1;              //因为周期边界条件， 如果穿越了（例如：右边出去，左边进去，那么计算活性取向应该-1*(r1-r4)）
                                        // ，在程序中写的活性速度方向要进行修正
        double labely = 1;               // 同x ，如果穿越边界，取向*-1
        double tempdr;                   //对取向归一化的 利用 dx/dr
        for (int i = 0; i < n; i++) {
            n3i = n3 + i;
            fCore = v0 / getNLco(dr, i);  // 活性速度=活性速度/周围粒子数
            //fCore = v0 ;
            labelx = Math.abs(dx[3 * n + i][i]) < L_2 ? 1 : -1;
            labely = Math.abs(dy[3 * n + i][i]) < L_2 ? 1 : -1;
            tempdr=Math.sqrt(dx[i][n3i]*dx[i][n3i]+dy[i][n3i]*dy[i][n3i]);  //计算dr=sqrt(dx^2 +dy^2)
            p[0] = labelx * dx[i][n3i] / tempdr;    //计算p_x=labelx * dx/dr
            p[1] = labely * dy[i][n3i] / tempdr;    //计算p_x=labelx * dx/dr
            insertRand(p);                          //对p进行旋转
            fActive[i][0] = fCore * p[0];           // 活性速度x方向分量=活性速度*p_x
            fActive[i][1] = fCore * p[1];           // 活性速度y方向分量=活性速度*p_y
        }
        return fActive;
    }

    public void insertRand(double[] p) {//加入自转 按照矩阵对应相乘
        final Random rnd = new Random();
        final double[] tmp = new double[2];
        double theata = Ts * rnd.nextGaussian();
        tmp[0] = Math.sin(theata) * p[1] + Math.cos(theata) * p[0];
        tmp[1] = -Math.sin(theata) * p[0] + Math.cos(theata) * p[1];
        p[0] = tmp[0];
        p[1] = tmp[1];
    }

    public double getNLco(double[][] dr, int i) {
        double N;
        final int n4 = n * 4;
        final double sigma = 1;
        final int n3 = (int) (3 * n);
        N = 0;
        for (int j = 0; j < n4; j++) {
            if ( cutr>=dr[i][j]) N++;

        }

        //if (N == 0) System.out.print(dr[i][i]);
        return N*4;
    }

    public double[] gexFy(double[][] fljcore, double[][] dy, int k) {
        double[] Fy;
        Fy = new double[n];
        final int n4 = n * 4;
        int nk = k * n;
        //int x,y;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n4; j++) {
                if(fljcore[i + nk][j]!=0)
                Fy[i] += fljcore[i + nk][j] * dy[i + nk][j];


            }
        }
        return Fy;

    }

    public double[] gexFx(double[][] fljcore, double[][] dx, int k) {
        double[] Fx;
        final int n4 = n * 4;
        int nk = k * n;
        Fx = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n4; j++) {
                if(0!=fljcore[i + nk][j])
                Fx[i] += fljcore[i + nk][j] * dx[i + nk][j];
            }
        }
        return Fx;
    }


    public double[][] getDr(double[][] dx, double[][] dy) {
        double[][] dr;
        double temp;
        dr = new double[4 * n][4 * n];
        final int n4 = n * 4;
        for (int i = 0; i < n4; i++) {
            for (int j = 0; j < i; j++) {
                //temp=Math.sqrt(dx[i][j] * dx[i][j] + dy[i][j]* dy[i][j]);
                dr[i][j] = (Math.abs(dx[i][j])<cutr)&&(Math.abs(dy[i][j])<cutr)?Math.sqrt(dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j]):999;
                //dr[i][j]=Math.sqrt(dx[i][j] * dx[i][j] + dy[i][j]* dy[i][j]);
            }
        }
        //

        for (int i = 0; i < n4; i++) {
            for (int j = i; j < n4; j++) {

                dr[i][j] = dr[j][i];
            }
        }


        return dr;
    }

    public double[][] getDy(double[][] Y) {
        double[][] dy;
        int temp1,temp2;
        final int n4 = n * 4;
        dy = new double[n4][n4];
        for (int i = 0; i < n4; i++) {
            temp1=i/n;
            temp2=i%n;
            for (int j = 0; j < i; j++) {

                dy[i][j] = Y[(temp1)][(temp2)] - Y[(j / n)][(j % n)];

            }

        }
        for (int i = 0; i < n4; i++) {
            for (int j = i; j < n4; j++) {
                dy[i][j] = -dy[j][i];

            }

        }
        return dy;
    }

    public double[][] getDx(double[][] X) {
        double[][] dx;
        final int n4 = n * 4;
        dx = new double[n4][n4];
        int temp1,temp2;
        for (int i = 0; i < n4; i++) {
            temp1=i/n;
            temp2=i%n;
            for (int j = 0; j < i; j++) {
                dx[i][j] = X[temp1][temp2] - X[j / n][j % n];
            }
        }
        for (int i = 0; i < n4; i++) {
            for (int j = i; j < n4; j++) {
                dx[i][j] = -dx[j][i];

            }

        }
        return dx;
    }

    public double[][] getFljcore(double[][] dr) {
        double[][] fcore;
        double ddr;
        final double e = 1;
        final double sigma = 1;
        final double sigma2 = 1/sigma/sigma;
        final int n4 = n * 4;
        double ddr6;
        final double cut = (24 * e * (2 * Math.pow(sigma / 0.8, 6) - 1) * Math.pow(sigma /  0.8, 8)) / sigma / sigma;
        fcore = new double[n4][n4];
        for (int i = 0; i < n4; i++) {
            for (int j = 0; j < i; j++) {




                if (dr[i][j] >1.112) {
                    fcore[i][j] = 0;
                    continue;
                }

                    //if ((i % n) == (j % n)) {
                        //fcore[i][j] = 0;
                       // continue;
                    //}//自己身上不受力

                    if (dr[i][j] > 0.8) {
                        ddr = sigma / dr[i][j];

                        fcore[i][j] = (24 * (2 * Math.pow(ddr, 6) - 1) * Math.pow(ddr, 8)) * sigma2;
                    } else {
                        fcore[i][j] = cut;
                    }


            }
        }


        for (int i = 0; i < n4; i++) {
            for (int j = i; j < n4; j++) {
                fcore[i][j] = fcore[j][i];
            }
        }


        return fcore;
    }

    public double[][][] getTcore(double[][] dr, double[][] dx, double[][] dy, double[][] X, double[][] Y) {
        double[][][] fcore;
        fcore = new double[4][n][2];
        double  attemp_ds=0.0000001;
        double temp1, temp2;
        double tempdx1, tempdy1, tempdr1;
        double tempdx2, tempdy2, tempdr2;
        double dx12, dy12, dx23, dy23, dx34, dy34;
        double E_bending;
        double softcore;
        for (int i = 0; i < n; i++) {
            tempdx1 = Math.abs(dx[n + i][i]) < L_2 ? dx[n + i][i] : Math.signum(X[0][i] - X[1][i]) * (L - Math.abs(X[1][i] - X[0][i]));
            tempdy1 = Math.abs(dy[n + i][i]) < L_2 ? dy[n + i][i] : Math.signum(Y[0][i] - Y[1][i]) * (L - Math.abs(Y[1][i] - Y[0][i]));
            tempdr1 = Math.sqrt(tempdy1 * tempdy1 + tempdx1 * tempdx1);
            temp1 = getT(tempdr1);
            fcore[0][i][0] = temp1 * (tempdx1);
            fcore[0][i][1] = temp1 * (tempdy1);
            dx12=tempdx1;
            dy12=tempdy1;
            tempdx2 = Math.abs(dx[2 * n + i][n + i]) < L_2 ? dx[2 * n + i][n + i] : Math.signum(X[1][i] - X[2][i]) * (L - Math.abs(X[2][i] - X[1][i]));
            tempdy2 = Math.abs(dy[2 * n + i][n + i]) < L_2 ? dy[2 * n + i][n + i] : Math.signum(Y[1][i] - Y[2][i]) * (L - Math.abs(Y[2][i] - Y[1][i]));
            tempdr2 = Math.sqrt(tempdy2 * tempdy2 + tempdx2 * tempdx2);
            temp2 = getT(tempdr2);
            dx23=tempdx2;
            dy23=tempdy2;

            //softcore=lam*(-tempdx2*tempdx1-tempdy1*tempdy2);
            //softcore=Math.abs(lam*(-tempdy2*tempdx1+tempdy1*tempdx2));

            fcore[1][i][0] = temp1 * (-tempdx1) + temp2 * (tempdx2);//dx[i][n+i] 即D12
            fcore[1][i][1] = temp1 * (-tempdy1) + temp2 * (tempdy2);//dy[i][n+i] 即D12
            tempdx1 = Math.abs(dx[3 * n + i][2 * n + i]) < L_2 ? dx[3 * n + i][2 * n + i] : Math.signum(X[2][i] - X[3][i]) * (L - Math.abs(X[3][i] - X[1][i]));
            tempdy1 = Math.abs(dy[3 * n + i][2 * n + i]) < L_2 ? dy[3 * n + i][2 * n + i] : Math.signum(Y[2][i] - Y[3][i]) * (L - Math.abs(Y[3][i] - Y[1][i]));
            tempdr1 = Math.sqrt(tempdy1 * tempdy1 + tempdx1 * tempdx1);
            temp1 = getT(tempdr1);
            dx34=tempdx1;
            dy34=tempdy1;

            //softcore=lam*(-tempdx2*tempdx1-tempdy1*tempdy2);

            fcore[2][i][0] = temp1 * (tempdx1) + temp2 * (-tempdx2);//dr[n+i][i] 即D23.d43
            fcore[2][i][1] = temp1 * (tempdy1) + temp2 * (-tempdy2);//dr[n+i][i] 即D21

            fcore[3][i][0] = temp1 * (-tempdx1);//dr[n+i][i] 即D21
            fcore[3][i][1] = temp1 * (-tempdy1);//dr[n+i][i] 即D21

            E_bending=getE(dx12, dy12, dx23, dy23, dx34, dy34);
            fcore[0][i][0]-=(getE(dx12+attemp_ds, dy12, dx23, dy23, dx34, dy34)-E_bending)/attemp_ds;
            fcore[0][i][1]-=(getE(dx12, dy12+attemp_ds, dx23, dy23, dx34, dy34)-E_bending)/attemp_ds;
            fcore[1][i][0]-=(getE(dx12-attemp_ds, dy12, dx23+attemp_ds, dy23, dx34, dy34)-E_bending)/attemp_ds;
            fcore[1][i][1]-=(getE(dx12, dy12-attemp_ds, dx23, dy23+attemp_ds, dx34, dy34)-E_bending)/attemp_ds;
            fcore[2][i][0]-=(getE(dx12, dy12, dx23-attemp_ds, dy23, dx34+attemp_ds, dy34)-E_bending)/attemp_ds;
            fcore[2][i][1]-=(getE(dx12, dy12, dx23, dy23-attemp_ds, dx34, dy34+attemp_ds)-E_bending)/attemp_ds;
            fcore[3][i][0]-=(getE(dx12, dy12, dx23, dy23, dx34-attemp_ds, dy34)-E_bending)/attemp_ds;
            fcore[3][i][1]-=(getE(dx12, dy12, dx23, dy23, dx34, dy34-attemp_ds)-E_bending)/attemp_ds;
        }
        return fcore;
    }

    public double getT(double r) {
        //final double k = 0.7;  //软链
        final double k = 5;
        final double r0 = 1.5;
        return k * (r - r0) / r;
    }
    public double getE(double dx12,double dy12,double dx23,double dy23,double dx34,double dy34) {
        return lam * ((dx12*dx23+dy12*dy23)/(Math.sqrt((dx12*dx12+dy12*dy12)*(dx23*dx23+dy23*dy23)))+(dx23*dx34+dy23*dy34)/(Math.sqrt((dx34*dx34+dy34*dy34)*(dx23*dx23+dy23*dy23))));
    }


}
