//2010/04/18
//rflab.

#ifndef _RF_EULER_ANGLES_TEMPLATE_H_
#define _RF_EULER_ANGLES_TEMPLATE_H_

#include <cassert>
#include <sstream>

#include "math.h"
#include "expressiontemplate.h"

//debug
#undef LOCALPRINTF
//#include <stdio.h>
//#define LOCALPRINTF(...) printf(__VA_ARGS__)
#define LOCALPRINTF(...)

///rock field!
namespace rf
{
	///���w
	namespace math
	{		
		///�I�C���[�p
		///pitch��heading��bank
		template<typename T>
		class TEulerAngles
		{
		public:
			T m_heading;
			T m_pitch;
			T m_bank;

		public:
			typedef T value_type;
			typedef TEulerAngles<T> self_type;
		public:

			///�f�X�g���N�^
			///virtual�͂��Ȃ�->���x���ٗl�ɗ�����
			~TEulerAngles(){}

			///�����Ȃ��R���X�g���N�^
			///���������܂߂ĉ������Ȃ�
			TEulerAngles()
			{
				LOCALPRINTF("TEulerAngles no arg  constructor\n");
			}

			///�A���O���w��R���X�g���N�^
			TEulerAngles(T heading_, T pitch_, T bank_)
			{
				LOCALPRINTF("TEulerAngles each member init constructor\n");
				m_heading = heading_;
				m_pitch   = pitch_;
				m_bank    = bank_;
			}


			///�R�s�[�R���X�g���N�^
			TEulerAngles(const self_type& rhs_) 
			{
				LOCALPRINTF("TEulerAngles copy constructor\n");
				m_heading = rhs_.m_heading;
				m_pitch = rhs_.m_pitch;
				m_bank = rhs_.m_bank;
			}

			///�R�s�[�R���X�g���N�^ for expression template
			template<class Lhs, class Op, class Rhs>
			TEulerAngles(const TEulerAnglesExpression<Lhs, Op, Rhs>& obj_)
			:m_heading(obj_.heading()),
			 m_pitch(obj_.pitch()),
			 m_bank(obj_.bank())
			{
				LOCALPRINTF("TEulerAngles expression copy constructor\n");
			}

			///�R�s�[
			self_type& operator = (const self_type& rhs_) 
			{
				LOCALPRINTF("TEulerAngles copy\n");
				m_heading = rhs_.m_heading;
				m_pitch   = rhs_.m_pitch;
				m_bank    = rhs_.m_bank;
				return *this;
			}				

			///�R�s�[ for expression template
			template<class Lhs, class Op, class Rhs>
			self_type& operator = (const TEulerAnglesExpression<Lhs, Op, Rhs>& rhs_) 
			{
				LOCALPRINTF("TEulerAngles expression copy\n");
				
				//A=A*B�̂悤�Ȏ���
				//�����̕ێ����Ă���l�𒼐ڂ�����Ɗ|���Z�̓r���ŗv�f�̒l������ւ���Ă��܂�����
				//�e���|�����I�u�W�F�N�g�Ɍv�Z���ʂ�ۑ����Ō�ɃX���b�v
				self_type obj(rhs_.heading(), rhs_.pitch(), rhs_.bank());
				*this = obj;
				return *this;
			}

			//-------------------------------------------------------subscript operators	
			value_type& heading(void)
			{
				return m_heading;
			}
				
			const value_type& heading(void) const
			{
				return m_heading;
			}
			
			value_type& pitch()
			{
				return m_pitch;
			}
		
			const value_type& pitch() const
			{
				return m_pitch;
			}
	
			value_type& bank()
			{
				return m_bank;
			}

			const value_type& bank() const
			{
				return m_bank;
			}


			////-------------------------------------------------------cast operators
			///�f�o�b�O�psting�L���X�g
			operator std::string () const
			{
				std::stringstream ss;
				ss << "euler angles\n\t";
				//ss << m_elem[0] << "+ (" << m_elem[1] << "*i + " << m_elem[2] << "*j + "<< m_elem[3] << "*k)";
				ss << "(" << m_heading << ", " <<m_pitch << ","<< m_bank << ")";
				return ss.str();
			}	

			//-------------------------------------------------------unary operators
			self_type operator - () 
			{
				return self_type(-m_heading, -m_pitch, -m_bank);
			}

			//-------------------------------------------------------compound assignment operators
			///�a
			self_type& operator += (const self_type& rhs_) 
			{
				m_heading += rhs_.m_heading;
				m_pitch   += rhs_.m_pitch;
				m_bank    += rhs_.m_bank;
				return *this;
			}

			///��
			self_type& operator -= (const self_type& rhs_) 
			{
				m_heading -= rhs_.m_heading;
				m_pitch   -= rhs_.m_pitch;
				m_bank    -= rhs_.m_bank;
				return *this;
			}

			///�X�J���{
			self_type& operator *= (const value_type& rhs_)
			{
				m_heading *= rhs_;
				m_pitch   *= rhs_;
				m_bank    *= rhs_;
				return *this;
			}
			
			//-------------------------------------------------------binary operators
			//�N���X�Œ�`����ꍇ�AExpressionTemplate�ɑΉ����l�����
			//����/�E��̋L�q�����ꂼ��ɕK�v�ŔώG�ɂȂ�
			//���̂��ߓ񍀉��Z�q�̓O���[�o���֐��Ƃ��Ē�`����

			//-------------------------------------------------------math operations
			///������
			///canonize:�i����j�`��F�߂�
			///���܂̂Ƃ���float only
			self_type& Canonize()
			{
				//pitch�𐳏��͈�[-��/2, ��/2]�ɕϊ�����
				//[-��, ��]�͈̔͂Ɏ��񂵂���ɗ��\�𔻒肵�Đ����͈͂ɕϊ�
				math::function::Cycle(&m_pitch, -definition::RF_PI, definition::RF_PI); 
				if (m_pitch < -definition::RF_OVER_2PI)
				{
					m_pitch   -= definition::RF_PI;
					m_heading += definition::RF_PI;
					m_bank    += definition::RF_PI;
				}
				else if(definition::OVER_PI2 < m_pitch)
				{
					m_pitch    = definition::RF_PI - m_pitch;
					m_heading += definition::RF_PI;
					m_bank    += definition::RF_PI;
				}

				//gimbal lock���Ă��Ȃ����
				//heading�𐳏��͈͂Ɏ���
				if (std::fabs(m_pitch)>definition::RF_OVER_2PI - definition::RF_EPSILON)
				{
					m_heading += m_bank;
					m_bank = 0.0f;
				}
				else
				{
					math::function::Cycle(&m_bank, -definition::RF_PI, definition::RF_PI); 
				}

				//heading�𐳏��͈͂Ɏ���
				math::function::Cycle(&m_heading, -definition::RF_PI, definition::RF_PI); 

				return *this;
			}
			
			///�����I�C���[�p�ɂȂ�
			self_type& Identity()
			{
				m_heading = 0;
				m_pitch   = 0;
				m_bank    = 0;
				return *this;
			}
		};
	}
}




#endif