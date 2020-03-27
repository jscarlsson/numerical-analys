# the first table

|             | N    | A\b                       | inv(A)*b                  |
|-------------|------|---------------------------|---------------------------|
| eiffel1.mat | 522  | 1.4606e-04 s = 146.06 �s  | 2.7411e-04 s = 274.41 �s  |
| eiffel2.mat | 798  | 4.6401e-04 s = 464.01 �s  | 0.0010 s = 1&nbsp;000 �s  |
| eiffel3.mat | 1122 | 0.0012 s = 1&nbsp;200 �s  | 0.0025 s = 2&nbsp;500 �s  |
| eiffel4.mat | 3184 | 0.0182 s = 18&nbsp;200 �s | 0.0444 s = 44&nbsp;400 �s |

# the second table

|             | Naiv         | LU          | Gles (ej LU) | Gles+LU    |
|-------------|--------------|-------------|--------------|------------|
| eiffel1.mat | 0.923002 s   | 0.304206 s  | 0.387095 s   | 0.037957 s |
| eiffel2.mat | 3.829454 s   | 0.902456 s  | 0.864064 s   | 0.073596 s |
| eiffel3.mat | 13.369861 s  | 2.235284 s  | 1.887173 s   | 0.168804 s |
| eiffel4.mat | 610.932927 s | 46.692336 s | 18.021481 s  | 2.161482 s |

*Varf�r g�r det snabbare att l�sa problemet med LU-faktorisering? - f�r det �r O(n<sup>3</sup>) ist�llet f�r O(n<sup>4</sup>)
*Vilken metod l�ser problemet snabbast? (Med/utan LU? Full/gles l�sare?) - LU+gles
*F�r vilken modell blir tidsvinsten st�rst? Varf�r? - den st�rsta, eftersom n<sup>3</sup> �kar mycket l�ngsammare �n n<sup>4</sup>