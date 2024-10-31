// vectorExpr.h
//

#ifndef VECTOREXPR_H
#define VECTOREXPR_H

#include "tag_traits.h"
#include <ostream>
#include <stdexcept>

namespace soundmath {
	struct vector_tag {};
	
	// operation tags
	template <class P>
	class Identity {
	public:
		static inline P apply (const P& m_a) {
			return m_a;
		}
	};
	template <class P>
	class UnaryMinus {
	public:
		static inline P apply (const P& m_a) {
			return -m_a;
		}
	};
	template <class P>
	class Exp {
	public:
		static inline P apply (const P& m_a) {
			return P (exp (m_a));
		}
	};
	#define DEFINE_EXPRESSION(ap,op) \
	template <class P> \
	class ap {\
	public:\
	  static inline P apply (const P& m_a, const P& m_b) { return m_a op m_b; }\
	};
	DEFINE_EXPRESSION (OpAdd, +)
	DEFINE_EXPRESSION (OpSub, -)
	DEFINE_EXPRESSION (OpMul, *)
	DEFINE_EXPRESSION (OpDiv, /)
	#undef DEFINE_EXPRESSION
	
	// expression tree nodes
	template <class P, class A, class F>
	class Expr1Func {
		A m_a;
	public:
		Expr1Func (const A& a) : m_a (a) {}
		P operator() (int n) const {
			return F::apply ( m_a (n) );
		}
		#ifdef SIZECHECK
			int size () const { return m_a.size (); }
		#endif		
	};
	template <class P, class A, class B, class Op>
	class Expr1BinOp {
		A m_a;
		B m_b;
	public:
		Expr1BinOp (const A& a, const B& b) : m_a (a), m_b (b) {}
		P operator() (int n) const {
			return Op::apply (m_a (n), m_b (n));
		}
		#ifdef SIZECHECK
			int size () const {
				if (m_a.size () != m_b.size ()) std::range_error (
					"size of Expr1BinOp operands should be the same");
				return m_a.size ();
			}
		#endif		
	};
	template <class P, class A, class Op>
	class Expr1OpScalar {
		A m_a;
		P m_b;
	public:
		Expr1OpScalar (const A& a, const P& b) : m_a (a), m_b (b) {}
		P operator() (int n) const {
			return Op::apply (m_a (n), m_b );
		}
		#ifdef SIZECHECK
			int size () const { return m_a.size (); }
		#endif			
	};
	template <class P, class B, class Op>
	class Expr1ScalarOp {
		P m_a;
		B m_b;
	public:
		Expr1ScalarOp (const P& a, const B& b) : m_a (a), m_b (b) {}
		P operator() (int n) const {
			return Op::apply (m_a, m_b (n) );
		}
		#ifdef SIZECHECK
			int size () const { return m_b.size (); }
		#endif		
	};
	template <class P, class V>
	class ConstRef1 {
		const V& m_v;
	public:
		ConstRef1 (const V& v) : m_v (v) {}
		P operator() (int n) const {
			return m_v (n);
		}
		#ifdef SIZECHECK
			int size () const { return m_v.size (); }
		#endif		
	};
	template <class P, class E>
	class Expr1 {
	private:
		E m_e;
	public:
		Expr1 (const E& e) : m_e (e) {}
		P operator() (int n) const {
			return m_e (n);
		}
		#ifdef SIZECHECK
			int size () const { return m_e.size (); }
		#endif		
	};
	
	// base class for vectors
	template <class P, class I>
	class VectorBase {
	public:
		typedef vector_tag tag_type;
		
		explicit VectorBase () {}
		int size () const {
			return static_cast<const I*> (this)->size ();
		}
		P operator() (int n) const {
			return static_cast<const I*> (this)->operator() (n);
		}
		template <class E> I& assignFrom (const Expr1<P, E>& x) {
			I *me = static_cast<I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size () ) std::range_error (
				"size should be the same for = operator");
		#endif
			for (int i = 0; i < me->size (); ++i) me->operator() (i) = x (i);
			return *me;
		}
		template <class V> I& assignFrom (const VectorBase<P, V>& x) {
			I *me = static_cast<I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size () ) std::range_error (
				"size should be the same for = operator");
		#endif
			for (int i = 0; i < me->size (); ++i) me->operator() (i) = x (i);
			return *me;
		}
		I& assignFrom (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->size (); ++i) me->operator() (i) = x;
			return *me;
		}
		I& assignFrom (const P* x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->size (); ++i) me->operator() (i) = x[i];
			return *me;
		}		
		template <class E> VectorBase<P, I>& operator+= (const Expr1<P, E>& x) {
			I* me = static_cast<I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size () ) std::range_error (
				"size should be the same for += operator");
		#endif
			for (int i = 0; i < me->size (); ++i) me->operator() (i) += x (i);
			return *me;
		}
		template <class V> VectorBase<P, I>& operator+= (const VectorBase<P, V>& x) {
			I *me = static_cast<I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size () ) std::range_error (
				"size should be the same for += operator");
		#endif
			for (int i = 0; i < me->size (); ++i) me->operator() (i) += x (i);
			return *me;
		}
		VectorBase<P, I>& operator+= (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->size (); ++i) me->operator() (i) += x;
			return *me;
		}
		template <class E> VectorBase<P,I >& operator-= (const Expr1<P, E>& x) {
			I* me = static_cast<I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size () ) std::range_error (
				"size should be the same for -= operator");
		#endif
			for (int i = 0; i < me->size (); ++i) me->operator() (i) -= x (i);
			return *me;
		}
		template <class V> VectorBase<P, I>& operator-= (const VectorBase<P, V>& x) {
			I *me = static_cast<I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size () ) std::range_error (
				"size should be the same for -= operator");
		#endif
			for (int i = 0; i < me->size (); ++i) me->operator() (i) -= x (i);
			return *me;
		}
		VectorBase<P, I>& operator-= (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->size (); ++i) me->operator() (i) -= x;
			return *me;
		}
		VectorBase<P, I>& operator*= (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->size (); ++i) me->operator() (i) *= x;
			return *me;
		}
		VectorBase<P, I>& operator/= (P x) {
			I *me = static_cast<I*> (this);
			for (int i = 0; i < me->size (); ++i) me->operator() (i) /= x;
			return *me;
		}
	
		template <class E> double inner (const Expr1<P, E>& x) const {
			double sum = 0.0;
			const I* me = static_cast<const I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size ()) std::range_error (
				"size should be the same for innner product");
		#endif
			for (int i = 0; i < me->size (); ++i) sum += me->operator() (i) * x (i);
			return sum;
		}
		template <class V> double inner (const VectorBase<P, V>& x) const {
			double sum = 0.0;
			const I* me = static_cast<const I*> (this);
		#ifdef SIZECHECK
			if (me->size () != x.size ()) std::range_error (
				"size should be the same for innner product");
		#endif
			for (int i = 0; i < me->size (); ++i) sum += me->operator() (i) * x (i);
			return sum;
		}
	};
	template <class T, class A>
	std::ostream& operator<< (std::ostream& s, const VectorBase<T, A>& m_a) {
		s << "1x" << m_a.size () << std::endl << "[";
		for (int i = 0; i< m_a.size (); ++i) {
			s << m_a (i);
			if (i != m_a.size () - 1) s << " ";
		}
	
		s << "]" << std::endl;
		return s;
	}
	
	// functions of VectorBase
	#define DEFINE_EXPRESSION(f,ap) \
	template <class P,class A> \
	Expr1<P, Expr1Func<P, ConstRef1<P,VectorBase<P, A> >, ap<P> > > \
	f(const VectorBase<P, A>& m_a) \
	{\
	   typedef Expr1Func<P, ConstRef1<P, VectorBase<P, A> >, ap<P> > ExprT;\
	   return Expr1<P, ExprT>(ExprT (ConstRef1<P,VectorBase<P, A> >(m_a))); \
	}
	DEFINE_EXPRESSION (ident, Identity)
	DEFINE_EXPRESSION (operator- , UnaryMinus)
	DEFINE_EXPRESSION (exp, Exp)
	#undef DEFINE_EXPRESSION
	
	// functions of Expr1
	#define DEFINE_EXPRESSION(f,ap) \
	template <class P, class A> \
	Expr1<P, Expr1Func<P, Expr1<P, A>, ap<P> > > \
	f(const Expr1<P, A>& m_a) \
	{\
	   typedef Expr1Func<P, Expr1<P, A> , ap<P> > ExprT;\
	   return Expr1<P, ExprT>(ExprT (m_a)); \
	}
	DEFINE_EXPRESSION (operator- , UnaryMinus)
	DEFINE_EXPRESSION (exp, Exp)
	#undef DEFINE_EXPRESSION
	
	// binary operations between Two VectorBases
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A,class B> \
	Expr1<P, Expr1BinOp<P, ConstRef1<P, VectorBase<P, A> >, ConstRef1<P,VectorBase<P, B> >, ap<P> > >\
	op (const VectorBase<P, A>& m_a, const VectorBase<P, B>& m_b) {\
	  typedef \
		Expr1BinOp<P, ConstRef1<P, VectorBase<P, A> >, ConstRef1<P, VectorBase<P, B> >, ap<P> > \
		  ExprT;\
	  return Expr1<P, ExprT>(ExprT (ConstRef1<P,VectorBase<P, A> >(m_a),\
					  ConstRef1<P,VectorBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	// binary operations between VectorBase and Scalar
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A>\
	Expr1<P, Expr1OpScalar<P, ConstRef1<P, VectorBase<P, A> >, ap<P> > > \
	op (const VectorBase<P, A>& m_a, const P& m_b) {\
	  typedef Expr1OpScalar<P, ConstRef1<P, VectorBase<P, A> >, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (ConstRef1<P,VectorBase<P, A> >(m_a),m_b));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	DEFINE_EXPRESSION (operator/, OpDiv)
	#undef DEFINE_EXPRESSION
	
	// binary operations between Scalar and VectorBase
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class B> \
	Expr1<P, Expr1ScalarOp<P, ConstRef1<P, VectorBase<P, B> >, ap<P> > > \
	op (const P& m_a, const VectorBase<P, B>& m_b) {\
	  typedef Expr1ScalarOp<P, ConstRef1<P, VectorBase<P, B> >, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a,ConstRef1<P,VectorBase<P, B> >(m_b)));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	#undef DEFINE_EXPRESSION
	
	// binary operations between Expr1 and Scalar
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A>\
	Expr1<P, Expr1OpScalar<P, Expr1<P, A>, ap<P> > > \
	op (const Expr1<P, A>& m_a, const P& m_b) {\
	  typedef Expr1OpScalar<P, Expr1<P, A>, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a,m_b));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	DEFINE_EXPRESSION (operator/, OpDiv)
	#undef DEFINE_EXPRESSION
	
	// binary operations between Scalar and Expr1
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class B> \
	Expr1<P, Expr1ScalarOp<P, Expr1<P, B>, ap<P> > > \
	op (const P& m_a, const Expr1<P, B>& m_b) {\
	  typedef Expr1ScalarOp<P, Expr1<P, B>, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a,m_b));\
	}
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	DEFINE_EXPRESSION (operator*, OpMul)
	#undef DEFINE_EXPRESSION
	
	// binary operations between VectorBase and Expr1
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A,class B>\
	Expr1<P, Expr1BinOp<P, ConstRef1<P, VectorBase<P, A> >, Expr1<P, B>, ap<P> > > \
	op (const VectorBase<P, A>& m_a, const Expr1<P, B>& m_b) {\
	  typedef Expr1BinOp<P, ConstRef1<P, VectorBase<P, A> >, Expr1<P, B>, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (ConstRef1<P,VectorBase<P, A> >(m_a), m_b));\
	}
	
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	// binary operations between Expr1 and VectorBase
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P,class A,class B>\
	Expr1<P, Expr1BinOp<P, Expr1<P, A>, ConstRef1<P, VectorBase<P, B> >, ap<P> > > \
	op (const Expr1<P, A>& m_a, const VectorBase<P, B>& m_b) {\
	  typedef Expr1BinOp<P, Expr1<P, A>, ConstRef1<P, VectorBase<P, B> >, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a,ConstRef1<P,VectorBase<P, B> >(m_b)));\
	}
	
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
	
	// binary operations between Two Expr1's
	#define DEFINE_EXPRESSION(op,ap) \
	template <class P, class A, class B>\
	Expr1<P, Expr1BinOp<P, Expr1<P, A>, Expr1<P, B>, ap<P> > >\
	op (const Expr1<P, A>& m_a, const Expr1<P, B>& m_b) {\
	  typedef Expr1BinOp<P, Expr1<P, A>, Expr1<P, B>, ap<P> > ExprT;\
	  return Expr1<P, ExprT>(ExprT (m_a, m_b));\
	}
	
	DEFINE_EXPRESSION (operator+, OpAdd)
	DEFINE_EXPRESSION (operator-, OpSub)
	#undef DEFINE_EXPRESSION
}

#endif // VECTOREXPR_H

// EOF


