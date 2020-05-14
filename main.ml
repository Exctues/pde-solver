(* Black Scholes solution for 1 dimensional case *)
open Owl;;

let e = 2.7182818284;;

let initial_condition is_call (s: float)  (k: float)  =
   if is_call then max (s -. k) 0.
    else max (k -. s) 0.;;

(* Useful function for debug *)
let print_matrix m = Mat.iteri (fun i j a -> Printf.printf "(%i,%i) %.1f\n" i j a) m;;

(*
 | creates a diagonal matrices from an array arr
 | and shifts by k depending on the sign:
 ^ k > 0 - to the upper right corner
 ^ k < 0 - to the lower left  corner
 ^ arr array which is used for diagonal matrix
 ^ n length of a given (Mat 1 n) array
*)
let diagonal arr n k =
   let a = Mat.zeros n n in
   for i = 0 to (n-1) do
     let value =  Mat.get a 1 i in
     let col = if k >= 0 then i else i-k in
     let row = if k >= 0 then i+k else i in
     Mat.set a col row  value;
   done;
   a;;

(*
  | returns C and D matrices for Crank-Nickolson scheme
  | In addition, returns coefficcients. Looks like C,D,a,b,c
  ^ sigma is volatility
  ^ m  is number pf points in price
  ^ r  is interest rate
  ^ dt is the timestep

*)
let construct_C_and_D sigma r dt m = 
  (* necessary coefficients *)
  let sig2 = sigma *. sigma in
  let a = Mat.zeros 1 (m+1) in
  let b = Mat.zeros 1 (m+1) in
  let c = Mat.zeros 1 (m+1) in

  for i = 0 to m do
    let fi = float_of_int i in
    let i2 = fi *. fi in
    Mat.set a 1 i (( dt   /.  4.)  *. (sig2 *. i2 -. r *. fi));
    Mat.set b 1 i ((-.dt  /.  2.)  *. (sig2 *. i2 +. r));
    Mat.set c 1 i (( dt   /.  4.)  *. (sig2 *. i2 +. r *. fi));
  done;
    
  (* Formulate necessary matrices *)
  let mat_C1 = Mat.neg (diagonal (Mat.get_slice_simple [[0];[2;m-1]] a) (m-2) (-1)) in
  let mat_C2_2_b = (Mat.get_slice_simple [[0];[1;m-1]] b) |> Mat.map (fun x -> 1. -. x) in
  let mat_C2 = (diagonal mat_C2_2_b (m - 1) 0) in 
  let mat_C3 = Mat.neg (diagonal (Mat.get_slice_simple [[0];[1;m-2]] c) (m-2) 1) in
  let mat_C = Mat.(mat_C1 + mat_C2 + mat_C3) in
  
  let mat_D1 = (diagonal (Mat.get_slice_simple [[0];[2;m-1]] a) (m-2) (-1)) in
  let mat_D2_2_b = (Mat.get_slice_simple [[0];[1;m-1]] b) |> Mat.map (fun x -> 1. +. x) in
  let mat_D2 = (diagonal mat_D2_2_b (m-1) 0) in
  let mat_D3 = (diagonal (Mat.get_slice_simple [[0];[1;m-2]] c) (m-2) 1) in
  let mat_D = Mat.(mat_D1 + mat_D2 + mat_D3) in
  mat_C, mat_D, a, b, c;;

let solve is_call r sigma matTime smin smax k n m =
  let ds = (smax -. smin) /. (float_of_int m) in
  let dt = matTime /. (float_of_int n) in
  let prices = Mat.zeros (m+2) (n+2) in

  (* set initial conditions *)
  for i = 1 to m + 1 do
    Mat.set prices i n (initial_condition is_call
      (float_of_int (i - 1) *. ds) k);
  done;
  
  (* Boundary condition *)
  for i = 1 to n + 1 do
    let common_factor = e ** (-.r *. (matTime -. (float_of_int (i-1)) *. dt)) in
    if is_call then
      begin 
        Mat.set prices 1 i 0.;
        Mat.set prices m i ((smax -. k) *. common_factor);
      end
    else
      begin
        Mat.set prices m i 0.;
        Mat.set prices 1 i ((k -. smin) *. common_factor);
      end
  done;

  (* Matrix preparation *)
  let mat_C, mat_D, a ,b ,c = construct_C_and_D sigma r dt m in

  let l, u, ipiv = Linalg.D.lu mat_C in

  let l_inv = Mat.inv l in
  let u_inv = Mat.inv u in

  let offset = Mat.zeros 1 m in
  
  (* Obtain a solution iteratively *)
  for i = (n-1) to 0 do
    Mat.set offset 1 0 (
      Mat.get a 0 2 *. (Mat.get prices 1 i +. Mat.get prices 1 (i+1)));
    Mat.set offset 1 0 (
      Mat.get c 0 (m-1) *. (Mat.get prices m i +. Mat.get prices m (i+1)));
    
    let px = Mat.get_slice_simple [[2;m];[i+1]] prices in 
    let x = Mat.((mat_D *@ px) + offset) in
    let x = Mat.(u_inv *@ (l_inv *@ x)) in

    Mat.set_slice_simple [[2;m];[i]] prices x;
  done;

  (* Show figures *)
  let useful_prices = Mat.get_slice_simple [[1;m+1];[1;n+1]] prices in
  let m_plus_one = m + 1 in
  let n_plus_one = n + 1 in
  let times = Mat.(sequential 1 (n_plus_one)
    |> map (fun x -> x -.1.)) in 
  let space_prices = Mat.(linspace 0. 10. (m_plus_one)
    |> map (fun x -> smin +. x *. ds)) in
  Plot.plot space_prices
   (Mat.get_slice_simple [[0;-1];[0]] useful_prices);
;;

let main = solve true 0.05 0.2 1. 0. 100. 60. 100 200;;