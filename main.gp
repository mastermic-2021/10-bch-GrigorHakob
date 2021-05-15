/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
x = ffgen(q,'a)
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);

syndrome(m, a, b) = {
	    alpha = 2*t -1;
	    return (sum(l=0,alpha , subst(m, 'x, a^(b+l)) * 'x^l));
	    }

 coder(m)={
	for(k = 0, q, 
	      a = ffprimroot(g);
	      for(b = 0, q-1,
	      	    info= List();
		    pade = bestapprPade(Mod(syndrome(msg,a,b),x^(2*t)));
		    D = denominator(pade);
		    N = numerator(pade);
		    
			  for(i = 0, q-2,
		    	   s= subst(D,'x,a^(-i));
			  if(s==0, val = subst((N/deriv(D))*(x^(b-1)), 'x , a^(-i));
			  t = fqx2int(val,g);
			  listput(info, t)));
	      
	      if(#info >5, return (Strchr(Vec(info))));
	      );
	);
}

recu = int2fqx(m, g);


secret=coder(recu);

print(secret);
