%X is a set of 8 3x1 vectors of normalized pixels that correspond to a set of 8 3x1 normalized vectors
%in Y. X is from the head image and Y is from the trailing image.
%Assumes: B'*E*A = 0.
function [R, t] = GetRandTfromE( E, An, Bn )

   [U, S, V] = svd(E);
   
   ss = [1,0,0;  %Put this back to the sigma1+sigma2/2 for dual identical singular values
        0,1,0;
        0,0,0];
  % ss = [ ( S(1,1) + S(2,2) ) / 2,                       0, 0;
  %                            0, ( S(1,1) + S(2,2) ) / 2, 0;
  %                            0,                       0, 0];    
  
   W = [ 0, -1,  0;
         1,  0,  0;
         0,  0,  1];

   tMat = U * W * ss * U';
   t1 = [ tMat(3,2), tMat(1,3), tMat(2,1) ]';
   t2 = -t1;
   
   R1 = U * W * V';
   R2 = U * inv(W) * V';
  
   
   
   Z = [  0,1,0;
         -1,0,0;
          0,0,0];      
   txx = U * Z * U';   
   tt1 = [ txx(3,2), txx(1,3), txx(2,1) ]';
   tt2 = -tt1;
   
   
   
   
   %There are 4 possible cases we need to test for:
   %Assume P = [I|0], and we need to find Pprime.
   
   %Case 1: Pp = [UWV' | +U3]
   yPrime = Bn(:,1);
   y = An(:,1);
   z3case1 = ( ( R2(1,:) - yPrime(1,1) * R2(3,:) ) * t1 ) / ...
             ( ( R2(1,:) - yPrime(1,1) * R2(3,:) ) * y );
      
   %Case 2: Pp = [UWV' | -U3]
   z3case2 = ( ( R2(1,:) - yPrime(1,1) * R2(3,:) ) * t2 ) / ...
             ( ( R2(1,:) - yPrime(1,1) * R2(3,:) ) * y );

   %Case 3: Pp = [UW'V' | +U3]
   z3case3 = ( ( R1(1,:) - yPrime(1,1) * R1(3,:) ) * t1 ) / ...
             ( ( R1(1,:) - yPrime(1,1) * R1(3,:) ) * y );
   
   %Case 4: Pp = [UW'V' | -U3]
   z3case4 = ( ( R1(1,:) - yPrime(1,1) * R1(3,:) ) * t2 ) / ...
             ( ( R1(1,:) - yPrime(1,1) * R1(3,:) ) * y );
            
   zcases = z3case1;
   zcases = horzcat( zcases, z3case2 );
   zcases = horzcat( zcases, z3case3 );
   zcases = horzcat( zcases, z3case4 );
   
   [ v, maxIdx ] = max( zcases );
   
   if( maxIdx == 1 )
       R = R2;
       t = t1;
   elseif( maxIdx == 2 )
       R = R2;
       t = t2;
   elseif( maxIdx == 3 )
       R = R1;
       t = t1
   else
       R = R1
       t = t2;
   end
end