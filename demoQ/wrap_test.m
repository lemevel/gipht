% compare two different ways of wrapping in C

% gcc wrap_test.c -lm -O3 -o wrap_test.maci64
% ./wrap_test.maci64 >! wrap_test.out

iar = load('wrap_test.out');
ii=iar(:,1);
aa=iar(:,2);
r1=iar(:,3);
r2=iar(:,4);
rd=iar(:,5);
r3=rwrapm2(aa);
figure;hold on;
plot(aa/2.0/pi,r1/2.0/pi,'r+-');
plot(aa/2.0/pi,r2/2.0/pi,'ko-');
plot(aa/2.0/pi,r3/2.0/pi,'bx-');
plot([min(aa/2.0/pi) max(aa/2.0/pi)], [+0.5  0.5], 'k--');
plot([min(aa/2.0/pi) max(aa/2.0/pi)], [-0.5 -0.5], 'k--');
%axis tight;
xlabel('unwrapped phase [cycles]');
ylabel('wrapped phase [cycles]');
legend('double in C','integer in C','rwrapm2 in Matlab');
printpdf(sprintf('%s.pdf',mfilename));
