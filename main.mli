(*
  | It solves BSM-PDE 1d case. 
  ^ r is interest rate,
  ^ sigma is volatility,
  ^ matTime is maturiy time (in years),
  ^ smin is min value of the underlying asset,
  ^ smax is max value of the underlying asset,
  ^ k is strike price of the underlying asset,
  ^ n is number of points in time,
  ^ m is number pf points in price.
*)

val solve : bool -> float -> float ->
  float -> float -> float -> float ->
  int -> int -> unit


(* 
  | return the value of initial condition border
  ^ k is strike price
  ^ s is price
*)
val initial_condition : bool -> float -> float -> float