import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class array {



    public static void main(String[] args)  {

        int m=1;
        int n=m*1;

        //double[] T={0.0,0.000005,0.00001,0.00002,0.00004,0.00008,0.00016,0.00032,0.00064,0.00128,0.00256,0.00512};
        double[]  lam=new double[m];
        double[]  Ts=new double[m];
        double[]  T=new double[m];
        double[]  v0=new double[m];
        for (int i =0;i<m;i++){T[i]=0.000001;Ts[i]=0;lam[i]=0*i;v0[i]=20;}
        //ExecutorService service=Executors.newFixedThreadPool(24);

        ThreadPoolExecutor service = (ThreadPoolExecutor) Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        service.setKeepAliveTime(10, TimeUnit.SECONDS);
        service.allowCoreThreadTimeOut(true);

        for(int i=0;i<n;i++){
            service.submit(new worm(i/m,T[i%m],Ts[i%m], lam[i%m],v0[i%m]));
            System.out.println("submit "+i);
        }
        System.out.println("submit finish");
        service.shutdown();

    }
}

