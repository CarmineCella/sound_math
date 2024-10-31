// matrixExpr.h
//

#ifndef MATRIXEXPR_H
#define MATRIXEXPR_H

#include "tag_traits.h"
#include "VecExpr.h"
#include <ostream>

namespace soundmath {
	struct matrix_tag {};
	
	template <class P, class A, class F>
	class Expr2Func {
		A m_a;
	public:
		Expr2Func (const A& a) : m_a (a) {}
		P operator() (int i, int j) const {
			return F::apply (m_a (i, j));
		}
	#ifdef SIZECHECK
		int rows () const {
			return m_a.rows ();
		}
		int cols () const {
			return m_a.cols ();
		}
	#endif
	};
	
	template <class P, class A, class B, class Op>
	class Expr2BinOp {
		A m_a;
		B m_b;
	public:
		Expr2BinOp (const A& a, const B& b) : m_a (a), m_b (b) {}
		P operator() (int i, int j) const {
			return Op::apply ( m_a (i, j), m_b (i, j));
		}
	};
	
	// Specialziation for Scalar Addition and Substractin
	// (Not Scalar Multiplication or Division)
	template <class P, class A, class Op>
	class Expr2OpScalar {
		A m_a;
		P m_b;
	public:
		Expr2OpScalar (const A& a, const P& b) : m_a (a), m_b (b) {}
		P operator() (int i, int j) const {
			return Op::apply (m_a (i, j), m_b);
		}
	};
	
	template <class P, class B, class Op>
	class Expr2ScalarOp {
		P m_a;
		B m_b;
	public:
		Expr2ScalarOp (const P& a, const B& b) : m_a (a), m_b (b) {}
		P operator() (int i, int j) const {
			return Op::apply (m_a, m_b (i, j));
		}
	};
	
	template <class P, class M>
	class ConstRef2 {
		const M& m_m;
	public:
		ConstRef2 (const M& m) : m_m (m) {}
		P operator() (int i, int j) const {
			return m_m (i, j);
		}
		int cols () const {
			return m_m.cols ();
		}
	};
	
	template <class P, class E>
	class Expr2 {
	private:
		E m_e;
	public:
		Expr2 (const E& e_) : m_e (e_) {}
		P operator() (int i, int j) const {
			return m_e (i, j);
		}
	};
	
	template <class P, class I>
	class MatrixBase {
	public:
		typedef matrix_tag tag_type;
		
		explicit MatrixBase() {}
		int size () const {
			return static_cast<const I*> (this)->size ();
		}
		int rows () const {
			return static_cast<const I*> (this)->rows ();
		}
		int cols () const {
			return static_cast<const I*> (this)->cols ();
		}
		P operator() (int i, int j) const {
			return static_cast<const I*> (this)->operator() (i, j);
		}
		template <class X> I& assignFrom (const Expr2<P, X>& rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) = rhs (i, j);
			}
			return *me;
		}
		template <class M> I& assignFrom (const MatrixBase<P, M>& x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) = x (i, j);
			}
			return *me;
		}
		I& assignFrom (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) = x;
			}
			return *me;
		}
		template <class X> MatrixBase<P, I>& operator+= (const Expr2<P, X>& rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) += rhs (i, j);
			}
			return *me;
		}
		template <class M> MatrixBase<P, I>& operator+= (const MatrixBase<P, M>& rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) += rhs (i, j);
			}
			return *me;
		}
		MatrixBase<P, I>& operator+= (P rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) += rhs;
			}
			return *me;
		}
		template <class X> MatrixBase<P, I>& operator-= (const Expr2<P, X>& rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) -= rhs (i, j);
			}
			return *me;
		}
		template <class M> MatrixBase<P, I>& operator-= (const MatrixBase<P, M>& rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) -= rhs (i, j);
			}
			return *me;
		}
		MatrixBase<P, I>& operator-= (P rhs) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) -= rhs;
			}
			return *me;
		}
		MatrixBase<P, I>& operator*= (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) *= x;
			}
			return *me;
		}
		MatrixBase<P, I>& operator/= (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->rows (); ++i) {
				for (int j = 0; j < me->cols (); ++j) me->operator() (i, j) /= x;
			}
			return *me;
		}
	};
	template <class T, class A>
	std::ostream& operator<< (std::ostream& s, const MatrixBase<T, A>& m_a) {
		s << m_a.rows () << "x" << m_a.cols () <<
			std::endl << "[" << std::endl;
		for (int i = 0; i< m_a.rows (); ++i) {
			s << "[";
			for (int j = 0; j < m_a.cols (); ++j) {
				s << m_a (i, j);
				if (j != m_a.cols () - 1) s << " ";
			}
	
			s << "]" << std::endl;
		}
		s << "]" << std::endl;
		return s;
	}
	
	// Matrix Vector multiplication
	// MatrixBase * VectorBase
	// MatrixBase * Expr1
	// Expr2 * VectorBase
	template <class P, class A, class B>
	class Expr1Reduct {
		A m_a;
		B m_b;
	public:
		Expr1Reduct (const A& a, const B& b) : m_a (a), m_b (b) {}
		P operator() (int i) const {
			P sum = m_a (i,0) * m_b (0);
			for (int j = 1; j < m_a.cols (); ++j) sum += m_a (i, j) * m_b (j);
			return sum;
		}
		#ifdef SIZECHECK
		int size () const {
			return m_b.size ();
		} 
		#endif
	};
	
	// Matrix Matrix Multiplication
	// MatrixBase * MatrixBase
	// MatrixBase * Expr2
	// Expr2 * MatrixBase
	// Expr2 * Expr2
	template <class P, class A, class B>
	class Expr2Reduct {
		A m_a;
		B m_b;
	public:
		Expr2Reduct (const A& a, const B& b) : m_a (a), m_b (b) {}
		P operator() (int i, int j) const {
			P sum = m_a (i, 0) * m_b (0, j);
			for (int k = 1; k < m_a.cols (); ++k) sum += m_a (i, k) * m_b (k, j);
			return sum;
		}
	};
	
	// Functions of MatrixBase
	#define DEFINE_EXPRESSION(f,ap) \
	template <class P, class A> \
	Expr2<P, Expr2Func<P, ConstRef2<P, MatrixBase<P, A> >, ap<P> > > \
	f(const MatrixBase<P, A>& m_a) \
	{\
	   typedef Expr2Func<P, ConstRef2<P, MatrixBase<P, A> >, ap<P> > ExprT;\
	   return Expr2<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a))); \
	}
	//DEFINE_EXPRESSION(ident, Identity)
	DEFINE_EXPRESSION (operator- , UnaryMinus)
	//DEFINE_EXPRESSION(exp, Exp)
	#undef DEFINE_EXPRESSION
	
	// Functions of Expr2
	#define DEFINE_EXPRESSION(f,ap) \
	template <class P, class E> \
	Expr2<P, Expr2Func<P, Expr2<P, E>, ap<P> > > \
	f(const Expr2<P, E>& m_a) \
	{\
	   typedef Expr2Func<P, Expr2<P, E>, ap<P> > ExprT;\
	   return Expr2<P, ExprT>(ExprT (m_a)); \
	}
	DEFINE_EXPRESSION (operator- , UnaryMinus)
	DEFINE_EXPRESSION (exp, Exp)
	#undef DEFINE_EXPRESSION
	
	//Binary operations between Two Dim2s
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P, class A, class B> \
	Expr2<P, Expr2BinOp<P, ConstRef2<P, MatrixBase<P, A> >, ConstRef2<P, MatrixBase<P, B> >, ap<P> > >\
	op (const MatrixBase<P, A>& m_a, const MatrixBase<P, B>& m_b) {\
	  typedef \
		Expr2BinOp<P, ConstRef2<P, MatrixBase<P, A> >, ConstRef2<P, MatrixBase<P, B> >, ap<P> > \
		  ExprT;\
	  return Expr2<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),\
					  ConstRef2<P, MatrixBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	// Multiplication with Two Dim2s
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr2<P, Expr2Reduct<P, ConstRef2<P, MatrixBase<P, A> >, ConstRef2<P, MatrixBase<P, B> > > >\
	op (const MatrixBase<P, A>& m_a, const MatrixBase<P, B>& m_b) {\
	  typedef \
		Expr2Reduct<P, ConstRef2<P, MatrixBase<P, A> >, ConstRef2<P, MatrixBase<P, B> > > \
		  ExprT;\
	  return Expr2<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),\
					  ConstRef2<P, MatrixBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION (operator*)
	#undef DEFINE_EXPRESSION
	
	//Binary operations between MatrixBase and Expr2
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A,class B> \
	Expr2<P, Expr2BinOp<P, ConstRef2<P, MatrixBase<P, A> >, Expr2<P, B>, ap<P> > >\
	op (const MatrixBase<P, A>& m_a, const Expr2<P, B>& m_b) {\
	  typedef Expr2BinOp<P, ConstRef2<P, MatrixBase<P, A> >, Expr2<P, B>, ap<P> > ExprT;\
	  return Expr2<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),m_b));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	/* Not Efficient !! // Muitiplication with  MatrixBase and Expr2
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr2<P, Expr2Reduct<P, ConstRef2<P, MatrixBase<P, A> >, Expr2<P, B> > >\
	op (const MatrixBase<P, A>& m_a, const Expr2<P, B>& m_b) {\
	  typedef \
		Expr2Reduct<P, ConstRef2<P, MatrixBase<P, A> >, Expr2<P, B> > \
		  ExprT;\
	  return Expr2<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),m_b));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	//Binary operations between MatrixBase and Scalar
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A> \
	Expr2<P, Expr2OpScalar<P, ConstRef2<P, MatrixBase<P, A> >, ap<P> > >\
	op (const MatrixBase<P, A>& m_a, P& m_b) {\
	  typedef \
		Expr2OpScalar<P, ConstRef2<P, MatrixBase<P, A> >, ap<P> > \
		  ExprT;\
	  return Expr2<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),m_b));\
	}
	// FIXME: are really used??
	DEFINE_EXPRESSION(operator+, OpAdd)
	DEFINE_EXPRESSION(operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	DEFINE_EXPRESSION (operator/, OpDiv)
	#undef DEFINE_EXPRESSION
	
	// impossible ! // Binary Operations between MatrixBase and VectorBase
	
	// Multiplication with MatrixBase and VectorBase
	#define DEFINE_EXPRESSION(op) \
	template <class P, class A, class B> \
	Expr1<P, Expr1Reduct<P, ConstRef2<P, MatrixBase<P, A> >, ConstRef1<P,VectorBase<P, B> > > >\
	op (const MatrixBase<P, A>& m_a, const VectorBase<P, B>& m_b) {\
	  typedef \
		Expr1Reduct<P, ConstRef2<P, MatrixBase<P, A> >, ConstRef1<P,VectorBase<P, B> > > \
		  ExprT;\
	  return Expr1<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),\
					  ConstRef1<P,VectorBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION (operator*)
	#undef DEFINE_EXPRESSION
	
	// impossible ! // Binary Operations between MatrixBase and Expr1
	
	/* Not Efficient !! // Muitiplicatino with  MatrixBase and Expr1
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr1<P, Expr1Reduct<P, ConstRef2<P, MatrixBase<P, A> >, Expr1<P, B> > >\
	op (const MatrixBase<P, A>& m_a, const Expr1<P, B>& m_b) {\
	  typedef \
		Expr1Reduct<P, ConstRef2<P, MatrixBase<P, A> >, Expr1<P, B> > \
		  ExprT;\
	  return Expr1<P, ExprT>(ExprT (ConstRef2<P, MatrixBase<P, A> >(m_a),m_b));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	/* Not Efficient! //Multiplication with Expr2 and MatrixBase
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr2<P, Expr2Reduct<P, Expr2<P, A>, ConstRef2<P, MatrixBase<P, B> > > >\
	op (const Expr2<P, A>& m_a, const MatrixBase<P, B>& m_b) {\
	  typedef \
		Expr2Reduct<P, Expr2<P, A>, ConstRef2<P, MatrixBase<P, B> > > \
		  ExprT;\
	  return Expr2<P, ExprT>(ExprT (m_a, ConstRef2<P, MatrixBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	//Binary operations between Two Expr2s
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P, class A, class B> \
	Expr2<P, Expr2BinOp<P, Expr2<P, A>, Expr2<P, B>, ap<P> > >\
	op (const Expr2<P, A>& m_a, const Expr2<P, B>& m_b) {\
	  typedef Expr2BinOp<P, Expr2<P, A>, Expr2<P, B>, ap<P> > ExprT;\
	  return Expr2<P, ExprT>(ExprT (Expr2<P, A>(m_a),Expr2<P, B>(m_b)));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	/* Not Efficient! //Multiplication with Two Expr2s
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P, class A, class B> \
	Expr2<P, Expr2Reduct<P, Expr2<P, A>, Expr2<P, B> > > >\
	op (const Expr2<P, A>& m_a, const Expr2<P, B>& m_b) {\
	  typedef Expr2Reduct<P, Expr2<P, A>, Expr2<P, B> > ExprT;\
	  return Expr2<P, ExprT>(ExprT (Expr2<P, A>(m_a),Expr2<P, B>(m_b)));\
	}
	DEFINE_EXPRESSION(operator+, OpAdd)
	DEFINE_EXPRESSION(operator-, OpSub)
	#undef DEFINE_EXPRESSION */
	
	//Binary operations between Expr2 and Scalar
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P, class A> \
	Expr2<P, Expr2OpScalar<P, Expr2<P, A>, ap<P> > >\
	op (const Expr2<P, A>& m_a, const P& m_b) {\
	  typedef Expr2OpScalar<P, Expr2<P, A>, ap<P> > ExprT;\
	  return Expr2<P, ExprT>(ExprT (m_a,m_b));\
	}
	// FIXME: are really used??
	DEFINE_EXPRESSION(operator+, OpAdd)
	DEFINE_EXPRESSION(operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	DEFINE_EXPRESSION (operator/, OpDiv)
	#undef DEFINE_EXPRESSION
	
	// impossible! //Binary Operations between Expr2 and VectorBase
	// impossible! //Multiplications with Expr2 and VectorBase
	
	// impossible! //Binary Operations between Expr2 and Expr1
	// impossible! //Multiplications with Expr2 and Expr1
	
	
	/* Not Efficient !! // Multiplication with Expr2 and VectorBase
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr1<P, Expr1Reduct<P, Expr2<P, A>, ConstRef1<P,VectorBase<P, B> > > >\
	op (const Expr2<P, A>& m_a, const VectorBase<P, B>& m_b) {\
	  typedef Expr1Reduct<P, Expr2<P, A>, ConstRef1<P,VectorBase<P, B> > > ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a,ConstRef1<P,VectorBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	//Binary operations between Expr2 and MatrixBase
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A,class B> \
	Expr2<P, Expr2BinOp<P, Expr2<P, A>, ConstRef2<P, MatrixBase<P, B> >, ap<P> > >\
	op (const Expr2<P, A>& m_a, const MatrixBase<P, B>& m_b) {\
	  typedef Expr2BinOp<P, Expr2<P, A>, ConstRef2<P, MatrixBase<P, B> >, ap<P> > ExprT;\
	  return Expr2<P, ExprT>(ExprT (m_a, ConstRef2<P, MatrixBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	/* Not Efficient! //Mutliplication with Expr2 and MatrixBase
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr2<P, Expr2Reduct<P, Expr2<P, A>, ConstRef2<P, MatrixBase<P, B> > > >\
	op (const Expr2<P, A>& m_a, const MatrixBase<P, B>& m_b) {\
	  typedef Expr2Reduct<P, Expr2<P, A>, ConstRef2<P, MatrixBase<P, B> > > ExprT;\
	  return Expr2<P, ExprT>(ExprT (m_a,ConstRef2<P, MatrixBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	// impossible ! // Binary Operations between Expr2 and VectorBase
	
	/* Not Efficient // Multiplication with Expr2 and VectorBase
	#define DEFINE_EXPRESSION(op) \
	template <class P, class A, class B> \
	Expr1<P, Expr1Reduct<P, Expr2<P, A>, ConstRef1<P,VectorBase<P, B> > > >\
	op (const Expr2<P, A>& m_a, const VectorBase<P, B>& m_b) {\
	  typedef \
		Expr1Reduct<P, Expr2<P, A>, ConstRef1<P,VectorBase<P, B> > > \
		  ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a, ConstRef1<P,VectorBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	// impossible ! // Binary Operations between Expr2 and Expr1
	
	/* Not Efficient !! // Muitiplicatino with  Expr2 and Expr1
	#define DEFINE_EXPRESSION(op) \
	template <class P,class A,class B> \
	Expr1<P, Expr1Reduct<P, Expr2<P, A>, Expr1<P, B> > >\
	op (const Expr2<P, A>& m_a, const Expr1<P, B>& m_b) {\
	  typedef \
		Expr1Reduct<P, Expr2<P, A>, Expr1<P, B> > \
		  ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a,m_b));\
	}
	DEFINE_EXPRESSION(operator*)
	#undef DEFINE_EXPRESSION */
	
	//Binary operations between Scalar and MatrixBase
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class B> \
	Expr2<P, Expr2ScalarOp<P, ConstRef2<P, MatrixBase<P, B> >, ap<P> > >\
	op (const P& m_a, const MatrixBase<P, B>& m_b) {\
	  typedef \
		Expr2ScalarOp<P, ConstRef2<P, MatrixBase<P, B> >, ap<P> > \
		  ExprT;\
	  return Expr2<P, ExprT>(ExprT (m_a,ConstRef2<P, MatrixBase<P, B> >(m_b)));\
	}
	// FIXME: are really used??
	DEFINE_EXPRESSION(operator+, OpAdd)
	DEFINE_EXPRESSION(operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	//DEFINE_EXPRESSION (operator/, OpDiv)	
	#undef DEFINE_EXPRESSION
	
	//Binary operations between Scalar and Expr2
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P, class B> \
	Expr2<P, Expr2ScalarOp<P, Expr2<P, B>, ap<P> > >\
	op (const P& m_a, const Expr2<P, B>& m_b) {\
	  typedef Expr2ScalarOp<P, Expr2<P, B>, ap<P> > ExprT;\
	  return Expr2<P, ExprT>(ExprT (m_a,m_b));\
	}
	// FIXME: are really used??
	DEFINE_EXPRESSION(operator+, OpAdd)
	DEFINE_EXPRESSION(operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	// DEFINE_EXPRESSION (operator/, OpDiv)	
	#undef DEFINE_EXPRESSION
	
	// ALREADY DEFINED! //Binary operations between Scalars
	// Already Defined! //Binary operations between Scalar and VectorBase
	// Already Defined! //Binary operations between Scalar and Expr1
	
	// impossible ! //Binary operations between VectorBase and MatrixBase
	// impossible ! //Multiplication with VectorBase and MatrixBase
	// impossible ! //Binary operations between VectorBase and Expr2
	// impossible ! //Multiplication with VectorBase and Expr2
	// Already Defined! //Binary operations between VectorBase and Scalar
	// Already Defined! //Binary operations between VectorBase and VectorBase
	// Already Defined! //Binary operations between VectorBase and Expr1
	
	// impossible ! //Binary operations between Expr1 and MatrixBase
	// impossible ! //Multiplication with Expr1 and MatrixBase
	// impossible ! //Binary operations between Expr1 and Expr2
	// impossible ! //Multiplication with Expr1 and Expr2
	// Already Defined! //Binary operations between Expr1 and Scalar
	// Already Defined! //Binary operations between Expr1 and VectorBase
	// Already Defined! //Binary operations between Expr1 and Expr1
}

#endif	// MATRIXEXPR_H

// EOF
