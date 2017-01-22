import java.awt.List;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;


public class SpecificGMM 
{	
	static ArrayList<Double> fileValues = new ArrayList<Double>(); 
	static double total = 0.0;
	static double count = 0.0;
	static double variance = 0.0;
	static double mean = 0.0;
	static double temp = 0.0;
	static int k = 3;
	static ArrayList<Double> uList = new ArrayList<Double>();
	static ArrayList<Double> prior;

	
//	public SpecificGMM(int k) {
//		SpecificGMM.k = k;
//	}


	public static void setPrior() {
		prior = new ArrayList<Double>(Collections.nCopies(k, 1.0/k));	
	}

	static ArrayList<Double> sigmaList = new ArrayList<Double>();
	
	double temp_u[] =new double[k];
	double temp_var[] = new double[k];
	static ArrayList<ArrayList<Double>> xgivenk = new ArrayList<ArrayList<Double>>(6000);
	static ArrayList<ArrayList<Double>> bgivenx = new ArrayList<ArrayList<Double>>(k);
	//static ArrayList<Double> sigmaList = new ArrayList<Double>();
	
	public static void setBgivenx()
	{
		ArrayList<Double> setzero = new ArrayList<Double>(Collections.nCopies(6000, 0.0));
		bgivenx.add(setzero);
	}

	
	public static void setXgivenk() {
		ArrayList<Double> setzero = new ArrayList<Double>(Collections.nCopies(6000, 0.0));
        xgivenk.add(setzero);
	}

	private static void readFile(String fileSource) throws IOException
	{
		 
		FileInputStream finput = new FileInputStream(fileSource);
		BufferedReader reader = new BufferedReader(new InputStreamReader(finput));	 
		String line = null;
		
		while ((line = reader.readLine()) != null)
		{			
			fileValues.add(Double.parseDouble(line));
		}	 
		//System.out.println(fileValues);
		reader.close();
	}
	
	
	
	public void setUSigmaList(ArrayList<Double> fileValues)
	{		
		Random r = new Random();
		double var = 1.0;
		//ArrayList<Double> priors = new ArrayList<>();
		for (int i = 0; i < k; i++) 
		{
			uList.add(fileValues.get(r.nextInt(fileValues.size())));
			sigmaList.add((i+1)*var);
		}	
		
	}
	
	public void xgivenk(double ulist,double sigList,int k)
	{
		
		ArrayList<Double> tempList = (ArrayList<Double>)xgivenk.get(k);
		int count=0;
		double temp1=0.0;
		double expTerm=0.0;
		double rootTerm=0.0;
		for(double val : fileValues)
		{	
			temp1 = -(val - ulist) * (val - ulist);
			temp1 = temp1/(2*sigList);
			expTerm = Math.exp(temp1);
			rootTerm = Math.sqrt(2*Math.PI*sigList);			
			tempList.set(count,expTerm/rootTerm);	
			count++;
		}
		xgivenk.set(k, tempList);
	
		
	}
	public void prior(int t)
	{
		ArrayList<Double> tempList = (ArrayList<Double>)bgivenx.get(t);
		double temp = 0;
		int i = 0;
		while(i < 6000)
		{
			temp = temp + tempList.get(i);
			i++;
		}
			prior.set(t, temp/6000);		
		
	}	
	
	public void bgivenx(int a)
	{
		ArrayList<Double> denom = new ArrayList<Double>(Collections.nCopies(6000, 0.0));
        ArrayList<Double> xkList = (ArrayList<Double>)xgivenk.get(a);
        ArrayList<Double> bkList = (ArrayList<Double>)bgivenx.get(a);
        setPrior();
		double den = 0,temp0 = 0;  
        
        for (int i = 0;i<k;i++)        	
        {
        	ArrayList<Double> tempList = xgivenk.get(i);
        	double p = prior.get(i);
        	for (int l = 0; l < 6000; l++) 
            {        
        		
        	den = denom.get(l);
        	double tmp1 = 	tempList.get(l);
        	
            temp0 = den + tmp1*p;
            denom.set(l, temp0);
            
            } 	
            
        }    

        int i = 0;		
        while(i<6000)
        {
        	double temp=0.0;
        	temp = bkList.get(i);
        	temp = temp + xkList.get(i)*prior.get(a)/denom.get(i);
        	bkList.set(i, temp);
        	i++;
        }
		bgivenx.set(a, bkList);
		
	}

	public void findUSigmaX(int m)
	{
	        double mTemp = 0, den = 0,sTemp = 0;
	        ArrayList<Double> meanlist = (ArrayList<Double>)bgivenx.get(m);
	        ArrayList<Double> sList = (ArrayList<Double>)bgivenx.get(m);
	        for(int j = 0;j<6000;j++){
	            mTemp=mTemp+meanlist.get(j)*fileValues.get(j);
	            den+=meanlist.get(j);
	        }        

	        uList.set(m, mTemp/den);
	       // System.out.println(uList);
	        for(int j = 0;j<6000;j++){
	            sTemp=sTemp+sList.get(j)*(fileValues.get(j)-uList.get(m))*(fileValues.get(j)-uList.get(m));
	            den+=sList.get(j);
	        }
	        sigmaList.set(m, sTemp/den);
	    }
	
	
	public static void main(String[] args) throws IOException
	{
		String fileSource;
		fileSource = args[0];
		//k = Integer.parseInt(args[1]);
		//int k = 2;
		
		SpecificGMM EM = new SpecificGMM();
		SpecificGMM.readFile(fileSource);
//		EM.setXgivenk();	
		EM.setUSigmaList(fileValues);	
		
		
		Double Mean = 0.0;
		int counter=0;
		boolean flag = true;
		
		while(flag)
		{
		for (int i = 0; i < k; i++) 
		{
			SpecificGMM.setXgivenk();
			EM.xgivenk(uList.get(i), sigmaList.get(i), i);
			
		}
		//System.out.println(xgivenk);
		for (int i = 0; i < k; i++) 
		{
			SpecificGMM.setBgivenx();
			EM.bgivenx(i);
		}
		
		for (int i = 0; i < k; i++) {
			EM.findUSigmaX(i);
		}
		
		for (int i = 0; i < k; i++) {
			EM.prior(i);
		}		
		xgivenk.clear();
		bgivenx.clear();		
		
		if(Mean/uList.get(0)>0.9999999 && Mean/uList.get(0)<1.0000001)
		{
			flag = false;
						
		}
		Mean = uList.get(0);
		counter++;	
		}
			
		
		System.out.println("Mean List :"+uList);
		System.out.println("Number of iterations :"+counter);
		System.out.println("Variance List :"+sigmaList);
		System.out.println("Number of iterations :"+counter);
        
		
		
	}
}