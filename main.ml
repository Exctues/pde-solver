(* Black Scholes solution for 1 dimensional case *)

(* Define Constants. *)
(* For convention all contants starts with _. *)

let _K  = 120 (* Strike price *)  
let _s0 = 100 (* Initial underlying asset price *)
let _Smax = 300 (* Max Value of the underlying asset *)
let _r  = 0.05 (* Risk free interest rate *)
let _sigma = 0.25 (* Volatility *)
let _NX = 10 (* number of space steps *)
let _NT = 10 (* number of time steps *)

(* We need to make *)
(* uniform partition of both space and time *)


(* Boundary and initial condition definition *)
let boundary_condition x y t = if x > y then 1.0 else 0.0
let initial_condition x y = 1.0

let sol  = Array.make_matrix _NX _NT 

(* Define a scheme to be used  *)

(* Implicit Euler *)




(* Visualization *)