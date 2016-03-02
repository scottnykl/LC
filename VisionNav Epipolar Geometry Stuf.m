%These are 8 matched points from the images (Scott's Matched Points)
%Ap =[
%   144,   253,   353,   139,   431, 325,   351,  248;
%    95,    91,   134,    25,   256,  96,   227,    2;
%    1 ,     1,     1,     1,     1,   1,     1,    1]
%                                                  
%                                                  
%Bp =[                                             
%   178,   252,   329,   167,   387, 309,   329,  248;
%    82,    72,   102,    25,   188,  75,   169,    5;
%     1,     1,     1,     1,     1,   1,     1,    1]

%Dan's Matched Points
Ap=[406 456 264 145 511 140 248 568;...
    101 248 236  96 85   27 141 14;...
     1   1   1   1   1   1    1  1];
                
Bp=[367 405 267 178 448 168 252 488;...
     77 183 177  83  61   27 109  8;...
      1   1   1   1   1   1   1   1];

ApNonMatch=[
   139, 140;
   209, 226;
     1,   1]

BpNonMatch=[
   559, 189;
    39, 114;
     1,   1]   
  
  
Klead = [
570.927050, 0.000000, 324.953512,
0.000000, 568.884811, 156.390491,
0., 0., 1.];

Ktrail = [
567.798343, 0.000000, 312.977357,
0.000000, 566.559477, 196.516183,
0., 0., 1.];

R = [
0.9996233464521829, 0.02734699290580352, -0.002303737560801904;
-0.02721660953720402, 0.9770668752974992, -0.2111861201918359;
-0.0035243996706521, 0.2111692761160844, 0.9774431520203978];

t = [-0.002; 0.062; 1.410];

[Ap Bp Klead Ktrail R t] = loadSceneConfiguration()

tx = [
    0, -t(3),  t(2);
 t(3),     0, -t(1);
-t(2),  t(1),     0]


E = tx * R;

F = inv( Ktrail') * E * inv(Klead);

%Part A-----------
EpConstraints = zeros(0);
EpConstraintsFund = zeros(0);
EpConstraintsFund_p = zeros(0);


%Pixels to Normalized
for m = 1:size(Ap,2)
   An(:,m) = inv(Klead) * Ap(:,m);
   Bn(:,m) = inv(Ktrail) * Bp(:,m);
end

%Check epipolar constraints
for m = 1:size(Ap,2)
   EpConstraints(:,m)= Bn(:,m)' * E * An(:,m);
end


%Use fundamental matrix instead of essential
for m = 1:size(Ap,2)
   EpConstraintsFund(:,m)= Bp(:,m)' * F * Ap(:,m);
end

EpConstraints
EpConstraintsFund

%Part B-----------

%Check epipolar constraints using the E matrix from the 8 point algorithm.
E_p = EightPtAlg( An, Bn );

for m = 1:size(Ap,2)
   EpConstraints8Pt(:,m)= Bn(:,m)' * E_p * An(:,m);
end

EpConstraints
EpConstraintsFund
EpConstraints8Pt

%Let's do a check and verify that the non matching points don't work
%well with either E or E_p

for m = 1:size(ApNonMatch,2)
   EpConstraintsNoMatch(:,m)= BpNonMatch(:,m)' * E * ApNonMatch(:,m);
end


for m = 1:size(ApNonMatch,2)
   EpConstraints8PtNoMatch(:,m)= BpNonMatch(:,m)' * E_p * ApNonMatch(:,m);
end

EpConstraintsNoMatch
EpConstraints8PtNoMatch

%Part C-----------
%[rr, tt] = GetRandTfromE( E )
%[rr, tt] = GetRandTfromE( E_p )


%For this to work we need a "good" E_p matrix.
%Create F_p matrix based on E_p matrix.
F_p = inv( Ktrail') * E_p * inv(Klead);
%Then convert normalized points to pixels to see error margins
for m = 1:size(Ap,2)
   EpConstraintsFund_p(:,m)= Bp(:,m)' * F_p * Ap(:,m);
end
%Now convert these to pixels

disp( 'Actual rpy')
DcmToRpy( R ) * 180/pi

disp( 'E rr,tt,rpy:')
[rr,tt] = GetRandTfromE( E, An, Bn )
DcmToRpy( rr ) * 180/pi
%WP = ProjectPoints( rr, tt, An, Bn )



disp( 'E_p rr,tt,rpy:')
[rr,tt] = GetRandTfromE( E_p, An, Bn )
DcmToRpy( rr ) * 180/pi
%WP = ProjectPoints( rr, tt, An, Bn )


