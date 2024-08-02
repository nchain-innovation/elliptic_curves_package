from typing import Optional
from copy import deepcopy

class BilinearPairing:
    def __init__(self, bilinear_pairing_curve, miller_output_type, easy_exponentiation, hard_exponentation):
        self.curve = bilinear_pairing_curve.curve
        self.twisted_curve = bilinear_pairing_curve.twisted_curve
        self.exp_t_minus_one = bilinear_pairing_curve.exp_t_minus_one
        self.miller_output_type = miller_output_type
        self.easy_exponentiation = easy_exponentiation
        self.hard_exponentiation = hard_exponentation

        return
    
    def miller_loop_on_base_curve(self, P, Q, denominator_elimination: Optional[str] = None):
        """
        Compute the Miller loop on P and Q on the base curve.
        This implementation is not optimised for the specific curve.

        The denominator_elimination variable decides which technique to use for denominator elimination
        """
        assert(denominator_elimination in [None, 'quadratic', 'cubic'])

        f = self.miller_output_type.identity()
        untwisted_Q = Q.to_base_curve()
        exp_t_minus_one = self.exp_t_minus_one

        if exp_t_minus_one[-1] == 1:
            T = untwisted_Q
        elif exp_t_minus_one[-1] == -1:
            T = -untwisted_Q
        else:
            raise ValueError('The most significant element of exp_t_minus_one must be non-zero')

        for i in range(len(exp_t_minus_one)-2,-1,-1):
            f = f.power(2)

            line_eval = T.line_evaluation(T,P)
            T = T + T

            match denominator_elimination:
                case 'quadratic':
                    line_eval = line_eval
                case 'cubic':
                    raise ValueError("To do!")
                case None:
                    line_eval = line_eval * T.line_evaluation(-T,P).invert()

            f = f.mul_by_line_eval(line_eval)
            
            if exp_t_minus_one[i] == 1:
                line_eval = T.line_evaluation(untwisted_Q,P)
                T = T + untwisted_Q
                
                match denominator_elimination:
                    case 'quadratic':
                        line_eval = line_eval
                    case 'cubic':
                        raise ValueError("To do!")
                    case None:
                        line_eval = line_eval * T.line_evaluation(-T,P).invert()
                
                f = f.mul_by_line_eval(line_eval)
            elif exp_t_minus_one[i] == -1:
                line_eval = T.line_evaluation(-untwisted_Q,P)
                T = T - untwisted_Q
                
                match denominator_elimination:
                    case 'quadratic':
                        line_eval = line_eval
                    case 'cubic':
                        raise ValueError("To do!")
                    case None:
                        line_eval = line_eval * T.line_evaluation(-T,P).invert()
                f = f.mul_by_line_eval(line_eval)
            else:
                pass

        return f
    
    def triple_miller_loop_on_base_curve(self, P1, P2, P3, Q1, Q2, Q3, denominator_elimination: Optional[str] = None):
        """
        Computes the product of three Miller loops on the base curve.
        It is not optimised at the moment (it computes three loops and multiplies them)
        """
        assert(denominator_elimination in [None, 'quadratic', 'cubic'])

        out1 = self.miller_loop_on_base_curve(P1,Q1,denominator_elimination)
        out2 = self.miller_loop_on_base_curve(P2,Q2,denominator_elimination)
        out3 = self.miller_loop_on_base_curve(P3,Q3,denominator_elimination)

        return out1 * out2 * out3
    
    def miller_loop_on_twisted_curve(self, P, Q, denominator_elimination: Optional[str] = None):
        """
        Compute the Miller loop on P and Q on the twisted curve.
        This implementation is not optimised for the specific curve.
        
        The denominator_elimination variable decides which technique to use for denominator elimination
        """
        assert(denominator_elimination in [None, 'quadratic', 'cubic'])

        f = self.miller_output_type.identity()
        twisted_P = P.to_twisted_curve()
        exp_t_minus_one = self.exp_t_minus_one

        if exp_t_minus_one[-1] == 1:
            T = deepcopy(Q)
        elif exp_t_minus_one[-1] == -1:
            T = -deepcopy(Q)
        else:
            raise ValueError('The most significant element of exp_t_minus_one must be non-zero')
        
        for i in range(len(exp_t_minus_one)-2,-1,-1):
            f = f.power(2)

            line_eval = T.line_evaluation(T,twisted_P)
            T = T + T

            match denominator_elimination:
                case 'quadratic':
                    line_eval = line_eval
                case 'cubic':
                    raise ValueError("To do!")
                case None:
                    line_eval = line_eval * T.line_evaluation(-T,twisted_P).invert()

            f = f.mul_by_line_eval(line_eval)
            
            if exp_t_minus_one[i] == 1:
                line_eval = T.line_evaluation(Q,twisted_P)
                T = T + Q
                
                match denominator_elimination:
                    case 'quadratic':
                        line_eval = line_eval
                    case 'cubic':
                        raise ValueError("To do!")
                    case None:
                        line_eval = line_eval * T.line_evaluation(-T,twisted_P).invert()
                
                f = f.mul_by_line_eval(line_eval)
            elif exp_t_minus_one[i] == -1:
                line_eval = T.line_evaluation(-Q,twisted_P)
                T = T - Q
                
                match denominator_elimination:
                    case 'quadratic':
                        line_eval = line_eval
                    case 'cubic':
                        raise ValueError("To do!")
                    case None:
                        line_eval = line_eval * T.line_evaluation(-T,twisted_P).invert()
                f = f.mul_by_line_eval(line_eval)
            else:
                pass

        return f
    
    def triple_miller_loop_on_twisted_curve(self, P1, P2, P3, Q1, Q2, Q3, denominator_elimination: Optional[str] = None):
        """
        Computes the product of three Miller loops on the base curve.
        It is not optimised at the moment (it computes three loops and multiplies them)
        """

        out1 = self.miller_loop_on_twisted_curve(P1,Q1,denominator_elimination)
        out2 = self.miller_loop_on_twisted_curve(P2,Q2,denominator_elimination)
        out3 = self.miller_loop_on_twisted_curve(P3,Q3,denominator_elimination)

        return out1 * out2 * out3
    
    def pairing(self, P, Q):
        """
        Computes the bilinear pairing on P and Q
        """
        if P.is_infinity() or Q.is_infinity():
            return self.miller_output_type.identity()
        else:
            out = self.miller_loop_on_base_curve(P,Q)
            out = self.easy_exponentiation(out)
            out = self.hard_exponentiation(out)

        return out

    def triple_pairing(self, P1, P2, P3, Q1, Q2, Q3):
        """
        Computes the product of three pairings

        The current implementation only allows the computation when neither of the Pi's or the Qi's is the point at infinity
        """

        assert(not(P1.is_infinity() or P2.is_infinity() or P3.is_infinity() or Q1.is_infinity() or Q2.is_infinity() or Q3.is_infinity())) 

        out = self.triple_miller_loop_on_base_curve(P1,P2,P3,Q1,Q2,Q3)
        out = self.easy_exponentiation(out)
        out = self.hard_exponentiation(out)

        return out
    


