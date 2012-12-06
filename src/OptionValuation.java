import org.apache.commons.math3.distribution.NormalDistribution;


public class OptionValuation {
	
	
	static double bankDaysPerYear = 252;
	static double rateDaysPerYear = 360;
	double T;
	
	
	double periods;
	
	double deltaT;
	
	
	double K;
	double S ;
	double SSVX3;
	double r = SSVX3;
	
	double sigma;
	
	double u;
	double d;
	
	double p;
	double p_1;
	
	double ePow_rt;
	
	public OptionValuation(double K, double S, double SSVX3, double sigma, int nrOfMonths){
		
		this.K = K;
		this.S = S;
		this.SSVX3 = SSVX3;
		this.r = SSVX3;
		this.sigma = sigma;
		this.T = bankDaysPerYear/12*nrOfMonths;
					
	}
	
	
	public void setBinConstants(int periods){
		this.periods = (double)periods;
		deltaT = T/periods/bankDaysPerYear;
		
		u = Math.exp(sigma*Math.sqrt(deltaT));
		d = 1/u;
		p = (Math.exp( (r*deltaT/periods)) - d ) / (u-d);
		p_1 = 1-p;
		ePow_rt = Math.exp(-deltaT*r);
	}
	
	public double callBlackAndScholeEuroOptionPricing(){
		double r = this.r/rateDaysPerYear;
		double T = this.T/bankDaysPerYear;
		
		double d1 = (Math.log(S/K) + (r + sigma*sigma/2)*(T))/(sigma*Math.sqrt(T));
		double d2 = d1 - sigma*Math.sqrt(T);
		
		
		NormalDistribution nd = new NormalDistribution();
		
		double C_st = nd.cumulativeProbability(d1)*S - nd.cumulativeProbability(d2)*K*Math.exp(-r*T);
		
		
		return C_st;
	}
	
	public double putBlackAndScholeEuroOptionPricing(){
		double r = this.r/rateDaysPerYear;
		double T = this.T/bankDaysPerYear;
	
		return K*Math.exp(-r*T) - S + callBlackAndScholeEuroOptionPricing();
	}
	
	
	
	
	public double callBinomialEuroOptionPricing(int periods){
		setBinConstants(periods);
		
		double vals[] = calculateFinalValuesCall();
		
		return calculateOptionPrice(vals);
	}
	
	public double putBinomialEuroOptionPricing(int periods){
		setBinConstants(periods);
		
		double vals[] = calculateFinalValuesPut();
		
		return calculateOptionPrice(vals);
	}
	
	
	
	
	public double callBinomialUSOptionPricing(int periods){
		return callBinomialEuroOptionPricing(periods);
	}
	
	public double putBinomialUSOptionPricing(int periods){
		setBinConstants(periods);
		
		double vals[] = calculateFinalValues();
		
		return calculateUSOptionPrice(vals);
	}
	
	
	
	private double[] calculateFinalValues(){
		double values[] = new double[(int)(periods+1)];
		
		for(int i = 0; i < periods+1; i++){
			values[i] = S*Math.pow(u, periods-i)*Math.pow(d, i);			
		}
		
		return values;
	}
	
	private double[] calculateFinalValuesPut(){
		double values[] = new double[(int)(periods+1)];
		
		for(int i = 0; i < periods+1; i++){
			double val = S*Math.pow(u, periods-i)*Math.pow(d, i);			
			
			if(K - val > 0){
				values[i] = K - val;
			}else{
				values[i] = 0;
			}
		}
		
		return values;
	}
	
	private double[] calculateFinalValuesCall(){
		double values[] = new double[(int)(periods+1)];
		
		for(int i = 0; i < periods+1; i++){
			double val = S*Math.pow(u, periods-i)*Math.pow(d, i);			
			
			if(val - K > 0){
				values[i] = val-K;
			}else{
				values[i] = 0;
			}
		}
		
		return values;
	}
	
	private double calculateOptionPrice(double values[]){
		double reduced[] = new double[values.length];
		System.arraycopy(values, 0, reduced, 0, values.length);

		int pointer = values.length;
		
		for(int i = 0; i < values.length+1; i++){
			
			for(int j = 0; j < pointer-1; j++){
				reduced[j] = ePow_rt*(p*reduced[j] + p_1*reduced[j+1]);
			}
			
			pointer--;
			
		}		
		return reduced[0];
	}
	
	private double calculateUSOptionPrice(double values[]){
		double usPut[] = new double[values.length];

		double usSell[] = new double[values.length];
		
		for(int i = 0; i < values.length; i++){
			double val = K - values[i];
			
			usSell[i] = val;
			
			if(val > 0){
				usPut[i] = val;

			}else{
				usPut[i] = 0;

			}
		}
		
		int pointer = values.length;
		
		for(int i = 0; i < values.length+1; i++){
			
			for(int j = 0; j < pointer-1; j++){
				
				values[j] = values[j]/u;
				
				usSell[j] = K - values[j]; 
				
				double alt1 = ePow_rt*(p*usPut[j] + p_1*usPut[j+1]);
				double alt2 = usSell[j];
				
				if(alt1 > alt2){
					usPut[j] = alt1;
				}else{
					usPut[j] = alt2;
				}
			}
			
			pointer--;
			
		}		
		
		
		return usPut[0];
	}
	
	
	public static void main(String args[]){
		System.out.println("Assignment 2-D\n");
		
		double BASprice = 6.3061112617;
		double K = 190;
		double S = 192.6;
		double SSVX3 = 0.0107;
		double sigma = 0.1203560368;
		int nrOfMonths = 3;
		int periods = 19;
		
	
//Uppgift 1
		System.out.println("---------------------------------------------------------------------\n");
		System.out.println("Part 1 data:");
		System.out.println("K:      " + K);
		System.out.println("S:      " + S);
		System.out.println("SSVX3:  " + SSVX3);
		System.out.println("sigma:  " + sigma);
		System.out.println("Months: " + nrOfMonths);
		System.out.println();
		
		
// Uppgift 1A		
		OptionValuation op = new OptionValuation(K, S, SSVX3, sigma, nrOfMonths);
		System.out.println("1A, Pricing for 3 month European call option Holmen: " + op.callBinomialEuroOptionPricing(10) );
		System.out.println();
		
// Uppgift 1B		
		double BinPrice = 0;
		BASprice = op.callBlackAndScholeEuroOptionPricing();
		
		int lastIndex = 0;
		for(int i = 1; i < 200; i++){
			BinPrice = op.callBinomialEuroOptionPricing(i);
			if( Math.abs((BASprice/BinPrice - 1)*100) > 0.99999){
				lastIndex = i;
			}
		}
		lastIndex++;
		
		BinPrice = op.callBinomialEuroOptionPricing(lastIndex);
		
		System.out.println("1B, Periods: " + lastIndex );
		
		System.out.println("    Binnomial price:      " + BinPrice);
		System.out.println("    Black & Schols price: " + BASprice);
		System.out.println("    diff:                 " + (BASprice/BinPrice - 1)*100 + "%\n");
		
		
//Uppgift 2
		
		BASprice = 6.3061112617;
		K = 170;
		S = 192.6;
		SSVX3 = 0.0107;
		sigma = 0.1203560368;
		nrOfMonths = 3;
		periods = 64;
		
		System.out.println("---------------------------------------------------------------------\n");
		System.out.println("Part 2 data:");
		System.out.println("K:      " + K);
		System.out.println("S:      " + S);
		System.out.println("SSVX3:  " + SSVX3);
		System.out.println("sigma:  " + sigma);
		System.out.println("months: " + nrOfMonths);		
		System.out.println();
// Uppgift 2A	
		
		
		
		op = new OptionValuation(K, S, SSVX3, sigma, nrOfMonths);
		double usPut = op.putBinomialUSOptionPricing(periods);
		System.out.println("2A, Pricing for 3 month American put option Holmen: " + usPut);
		System.out.println();
		
//Uppgift 2B		
		
		double euPut = op.putBinomialEuroOptionPricing(periods);

		System.out.println("2A, Diffrence of American and Europeen Put: ");
		System.out.println("    American Put:  " + usPut);
		System.out.println("    European Put: -" + euPut);
		System.out.println("    Result:        " + (usPut-euPut));
		

//Uppgift 3
		
		

	}
}
