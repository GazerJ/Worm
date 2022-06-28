import java.io.*;


public class np {

    public static double[]  loadtxt(String path,int n) throws IOException {
        File file =new File(path);
        double[] data=new double[n];
        try {
            FileReader fis=new FileReader(file);
            BufferedReader dis=new BufferedReader(fis);
            String tmp=null;
            int i=0;
            while((tmp= dis.readLine())!=null){
                data[i]=Double.parseDouble(tmp);
                i++;
                if(i>=n)break;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return  data;
    }



    public static void  savetxt(String path,double[] data) throws IOException {
        File file =new File(path);
        try {
            FileWriter fis=new FileWriter(file);
            BufferedWriter dis=new BufferedWriter(fis);
            for(int i=0;i< data.length;i++) {
                dis.write(Double.toString(data[i]));
                dis.newLine();
            }
            dis.close();
            fis.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void  savetxt(String path,double[][] data) throws IOException {
        File file =new File(path);
        try {
            FileWriter fis=new FileWriter(file);
            BufferedWriter dis=new BufferedWriter(fis);
            for(int i=0;i< data.length;i++) {
                for(int j=0;j<data[i].length;j++) {
                    dis.write(Double.toString(data[i][j]));
                    dis.write("\t");
                }
                dis.newLine();
            }
            dis.close();
            fis.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void  savetxtN(String path,double[][] data,int n) throws IOException {
        File file =new File(path);
        try {
            FileWriter fis=new FileWriter(file);
            BufferedWriter dis=new BufferedWriter(fis);
            for(int i=0;i< data.length/n;i++) {
                for(int j=0;j<data[i].length;j++) {
                    dis.write(Double.toString(data[i*n][j]));
                    dis.write("\t");
                }
                dis.newLine();
            }
            dis.close();
            fis.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void savetxtCutAdd(String filepath , double[][] data,int n) {
        FileWriter fw = null;
        try {
            //如果文件存在，则追加内容；如果文件不存在，则创建文件
            File f=new File(filepath);
            fw = new FileWriter(f, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
        PrintWriter pw = new PrintWriter(fw);

        for(int i=0;i< data.length/n;i++) {
            for(int j=0;j<data[i].length;j++) {
                pw.print(Double.toString(data[i*n][j]));
                pw.print("\t");
            }
            pw.print("\n");
        }

        pw.flush();
        try {
            fw.flush();
            pw.close();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }





    public static void fill(double[][] x1, double[][] x2, double[][] x3, double[][] x4, double[][] y1, double[][] y2, double[][] y3, double[][] y4) {
    for(int i=0;i< x1.length;i++){
        for(int j=0;j< x1[0].length;j++){

            x1[i][j]=0;
            x2[i][j]=0;
            x3[i][j]=0;
            x4[i][j]=0;
            y1[i][j]=0;
            y2[i][j]=0;
            y3[i][j]=0;
            y4[i][j]=0;


        }

    }

    }
}
