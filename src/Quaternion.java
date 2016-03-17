import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/** Quaternions. Basic operations. */
public class Quaternion {

	// TODO!!! Your fields here!

	/**
	 * Constructor from four double values.
	 * 
	 * @param a
	 *            real part
	 * @param b
	 *            imaginary part i
	 * @param c
	 *            imaginary part j
	 * @param d
	 *            imaginary part k
	 */

	private double ldr, ldi, ldj, ldk;

	//private int hashCode;

	// used to compare to floats, if difference is less than DELTA, floats are
	// considered to be equal
	public static final double DELTA = 0.000001;

	public Quaternion(double a, double b, double c, double d) {
		// TODO!!! Your constructor here!
		this.ldr = a;
		this.ldi = b;
		this.ldj = c;
		this.ldk = d;
	}

	/**
	 * Real part of the quaternion.
	 * 
	 * @return real part
	 */
	public double getRpart() {
		return this.ldr; // TODO!!!
	}

	/**
	 * Imaginary part i of the quaternion.
	 * 
	 * @return imaginary part i
	 */
	public double getIpart() {
		return this.ldi; // TODO!!!
	}

	/**
	 * Imaginary part j of the quaternion.
	 * 
	 * @return imaginary part j
	 */
	public double getJpart() {
		return this.ldj; // TODO!!!
	}

	/**
	 * Imaginary part k of the quaternion.
	 * 
	 * @return imaginary part k
	 */
	public double getKpart() {
		return this.ldk; // TODO!!!
	}

	/**
	 * Conversion of the quaternion to the string.
	 * 
	 * @return a string form of this quaternion: "a+bi+cj+dk" (without any
	 *         brackets)
	 */
	@Override
	public String toString() {
		return String.format("%.0f%+.0fi%+.0fj%+.0fk", this.ldr, this.ldi, this.ldj, this.ldk); // TODO!!!
	}

	/**
	 * Conversion from the string to the quaternion. Reverse to
	 * <code>toString</code> method.
	 * 
	 * @throws IllegalArgumentException
	 *             if string s does not represent a quaternion (defined by the
	 *             <code>toString</code> method)
	 * @param s
	 *            string of form produced by the <code>toString</code> method
	 * @return a quaternion represented by string s
	 */
	public static Quaternion valueOf(String s) {
		double a, b, c, d;

		String pattern = "([-+]{0,1}(\\d)+)+([-+]{1,1}(\\d)+)i+([-+]{1,1}(\\d)+)j+([-+]{1,1}(\\d)+)k+";
		Pattern pat = Pattern.compile(pattern);
		Matcher mat = pat.matcher(s);
		if (!mat.matches() || mat.groupCount() != 8) {
			throw new IllegalArgumentException(String.format("%s does not conform to Quaternion string format 'a+bi+cj+dk'", s));

		}

		try {
			a = Double.valueOf(mat.group(1));
			b = Double.valueOf(mat.group(3));
			c = Double.valueOf(mat.group(5));
			d = Double.valueOf(mat.group(7));

		} catch (RuntimeException e) {
			throw new IllegalArgumentException(
					String.format("%s does not conform to Quaternion string format 'a+bi+cj+dk' (%s)", s, e.getMessage()));
		}

		return new Quaternion(a, b, c, d); // TODO!!!

	}

	/**
	 * Clone of the quaternion.
	 * 
	 * @return independent clone of <code>this</code>
	 */
	@Override
	public Object clone() throws CloneNotSupportedException {

		return new Quaternion(this.getRpart(), this.getIpart(), this.getJpart(), this.getKpart()); // TODO!!!
	}

	/**
	 * Test whether the quaternion is zero.
	 * 
	 * @return true, if the real part and all the imaginary parts are (close to)
	 *         zero
	 */
	public boolean isZero() {
		return (Math.abs(getRpart()) < DELTA && 
				Math.abs(getIpart()) < DELTA && 
				Math.abs(getJpart()) < DELTA && 
				Math.abs(getKpart()) < DELTA); // TODO!!!
	}

	/**
	 * Conjugate of the quaternion. Expressed by the formula
	 * conjugate(a+bi+cj+dk) = a-bi-cj-dk
	 * 
	 * @return conjugate of <code>this</code>
	 */
	public Quaternion conjugate() {
		return new Quaternion(
			getRpart()
			, -getIpart()
			, -getJpart()
			, -getKpart()
				);
	}

	/**
	 * Opposite of the quaternion. Expressed by the formula opposite(a+bi+cj+dk)
	 * = -a-bi-cj-dk
	 * 
	 * @return quaternion <code>-this</code>
	 */
	public Quaternion opposite() {
		
		return new Quaternion(
				-getRpart()
				, -getIpart()
				, -getJpart()
				, -getKpart()
				); // TODO!!!
	}

	/**
	 * Sum of quaternions. Expressed by the formula (a1+b1i+c1j+d1k) +
	 * (a2+b2i+c2j+d2k) = (a1+a2) + (b1+b2)i + (c1+c2)j + (d1+d2)k
	 * 
	 * @param q
	 *            addend
	 * @return quaternion <code>this+q</code>
	 */
	public Quaternion plus(Quaternion q) {
		
		return new Quaternion(
				this.getRpart()+q.getRpart()
				,this.getIpart()+q.getIpart()
				,this.getJpart()+q.getJpart()
				,this.getKpart()+q.getKpart()
				)	; // TODO!!!
	}

	/**
	 * Product of quaternions. Expressed by the formula (a1+b1i+c1j+d1k) *
	 * (a2+b2i+c2j+d2k) = (a1a2-b1b2-c1c2-d1d2) + (a1b2+b1a2+c1d2-d1c2)i +
	 * (a1c2-b1d2+c1a2+d1b2)j + (a1d2+b1c2-c1b2+d1a2)k
	 * 
	 * @param q
	 *            factor
	 * @return quaternion <code>this*q</code>
	 */
	public Quaternion times(Quaternion q) {
		return new Quaternion(
					getRpart()*q.getRpart() - getIpart()*q.getIpart() - getJpart()*q.getJpart() - getKpart()*q.getKpart()
				,   getRpart()*q.getIpart() + getIpart()*q.getRpart() + getJpart()*q.getKpart() - getKpart()*q.getJpart()
				,   getRpart()*q.getJpart() - getIpart()*q.getKpart() + getJpart()*q.getRpart() + getKpart()*q.getIpart()
				,   getRpart()*q.getKpart() + getIpart()*q.getJpart() - getJpart()*q.getIpart() + getKpart()*q.getRpart()
				);
	}

	/**
	 * Multiplication by a coefficient.
	 * 
	 * @param r
	 *            coefficient
	 * @return quaternion <code>this*r</code>
	 */
	public Quaternion times(double r) {
		return new Quaternion(
				r*getRpart()
				, r*getIpart()
				, r*getJpart()
				, r*getKpart()
				);
	}

	/**
	 * Inverse of the quaternion. Expressed by the formula 1/(a+bi+cj+dk) =
	 * a/(a*a+b*b+c*c+d*d) + ((-b)/(a*a+b*b+c*c+d*d))i +
	 * ((-c)/(a*a+b*b+c*c+d*d))j + ((-d)/(a*a+b*b+c*c+d*d))k
	 * 
	 * @return quaternion <code>1/this</code>
	 */
	public Quaternion inverse() {
		
		if (isZero()) {
			throw new RuntimeException("Division by zero!");
		}
		
		double d = getRpart()*getRpart() +
				getIpart()*getIpart() + 
				getJpart()*getJpart() +
				getKpart()*getKpart();
		return new Quaternion(
				(getRpart() / d )
				, ((-getIpart()) / d )
				, ((-getJpart()) / d )
				, ((-getKpart()) / d )
				);
				
	}

	/**
	 * Difference of quaternions. Expressed as addition to the opposite.
	 * 
	 * @param q
	 *            subtrahend
	 * @return quaternion <code>this-q</code>
	 */
	public Quaternion minus(Quaternion q) {
		return plus(q.opposite());
	}

	/**
	 * Right quotient of quaternions. Expressed as multiplication to the
	 * inverse.
	 * 
	 * @param q
	 *            (right) divisor
	 * @return quaternion <code>this*inverse(q)</code>
	 */
	public Quaternion divideByRight(Quaternion q) {
		
		if (q.isZero()) {
			throw new RuntimeException("Division by zero!");
		}
		
		return times(q.inverse());

				
	}

	/**
	 * Left quotient of quaternions.
	 * 
	 * @param q
	 *            (left) divisor
	 * @return quaternion <code>inverse(q)*this</code>
	 */
	public Quaternion divideByLeft(Quaternion q) {
		
		if (q.isZero()) {
			throw new RuntimeException("Division by zero!");
		}
		
		return q.inverse().times(this);
	}

	private boolean doubleEquals(double d1, double d2) {
		return Math.abs(d1 - d2) < DELTA;
	}
	
	/**
	 * Equality test of quaternions. Difference of equal numbers is (close to)
	 * zero.
	 * 
	 * @param qo
	 *            second quaternion
	 * @return logical value of the expression <code>this.equals(qo)</code>
	 */
	@Override
	public boolean equals(Object qo) {
		if (qo == null)
			return false;
		if (qo.getClass() != this.getClass())
			return false;

		Quaternion oo = (Quaternion) qo;
		return (
				doubleEquals(getRpart(), oo.getRpart()) &&
				doubleEquals(getIpart(), oo.getIpart()) &&
				doubleEquals(getJpart(), oo.getJpart()) &&
				doubleEquals(getKpart(), oo.getKpart()) 
				
				);

	}

	/**
	 * Dot product of quaternions. (p*conjugate(q) + q*conjugate(p))/2
	 * 
	 * @param q
	 *            factor
	 * @return dot product of this and q
	 */
	public Quaternion dotMult(Quaternion q) {
		
		return this.times(q.conjugate()).plus(q.times(this.conjugate())).times(0.5);
	}

	/**
	 * Integer hashCode has to be the same for equal objects.
	 * 
	 * @author NASA WorldWind
	 * @see <a href="http://worldwind31.arc.nasa.gov/svn/trunk/WorldWind/src/gov/nasa/worldwind/geom/Quaternion.java">http://worldwind31.arc.nasa.gov/svn/trunk/WorldWind/src/gov/nasa/worldwind/geom/Quaternion.java</a>
	 * 
	 * @return hashcode
	 */
	@Override
	public int hashCode() {
		
            int result;
            long tmp;
            tmp = Double.doubleToLongBits(this.ldr);
            result = (int) (tmp ^ (tmp >>> 32));
            tmp = Double.doubleToLongBits(this.ldi);
            result = 31 * result + (int) (tmp ^ (tmp >>> 32));
            tmp = Double.doubleToLongBits(this.ldj);
            result = 31 * result + (int) (tmp ^ (tmp >>> 32));
            tmp = Double.doubleToLongBits(this.ldk);
            result = 31 * result + (int) (tmp ^ (tmp >>> 32));
        
        return result;
	}

	/**
	 * Norm of the quaternion. Expressed by the formula norm(a+bi+cj+dk) =
	 * Math.sqrt(a*a+b*b+c*c+d*d)
	 * 
	 * @return norm of <code>this</code> (norm is a real number)
	 */
	public double norm() {
		double d = getRpart()*getRpart() +
				getIpart()*getIpart() + 
				getJpart()*getJpart() +
				getKpart()*getKpart();
		return Math.sqrt(d);
	}

	/**
	 * Main method for testing purposes.
	 * 
	 * @param arg
	 *            command line parameters
	 */
	public static void main(String[] arg) {

		Quaternion arv1 = new Quaternion(-1., 1, 2., -2.);
		if (arg.length > 0)
			arv1 = valueOf(arg[0]);
		System.out.println("first: " + arv1.toString());
		System.out.println("real: " + arv1.getRpart());
		System.out.println("imagi: " + arv1.getIpart());
		System.out.println("imagj: " + arv1.getJpart());
		System.out.println("imagk: " + arv1.getKpart());
		System.out.println("isZero: " + arv1.isZero());
		System.out.println("conjugate: " + arv1.conjugate());
		System.out.println("opposite: " + arv1.opposite());
		System.out.println("hashCode: " + arv1.hashCode());

		Quaternion k1 = new Quaternion(3., 7., 5., 2.);
		Quaternion k2 = new Quaternion(3., -9., 5., 2.);
		Quaternion k3 = new Quaternion(3., 7., 5., 2.);

		System.out.println(k1);
		System.out.println(k2);
		System.out.println(k3);

		if (k1.equals(k2)) {
			System.out.println("EQUALS");
		}

		// Quaternion res = null;
		// try {
		// res = (Quaternion)arv1.clone();
		// } catch (CloneNotSupportedException e) {};
		// System.out.println ("clone equals to original: " + res.equals
		// (arv1));
		// System.out.println ("clone is not the same object: " + (res!=arv1));
		// System.out.println ("hashCode: " + res.hashCode());
		// res = valueOf (arv1.toString());
		// System.out.println ("string conversion equals to original: "
		// + res.equals (arv1));
		// Quaternion arv2 = new Quaternion (1., -2., -1., 2.);
		// if (arg.length > 1)
		// arv2 = valueOf (arg[1]);
		// System.out.println ("second: " + arv2.toString());
		// System.out.println ("hashCode: " + arv2.hashCode());
		// System.out.println ("equals: " + arv1.equals (arv2));
		// res = arv1.plus (arv2);
		// System.out.println ("plus: " + res);
		// System.out.println ("times: " + arv1.times (arv2));
		// System.out.println ("minus: " + arv1.minus (arv2));
		// double mm = arv1.norm();
		// System.out.println ("norm: " + mm);
		// System.out.println ("inverse: " + arv1.inverse());
		// System.out.println ("divideByRight: " + arv1.divideByRight (arv2));
		// System.out.println ("divideByLeft: " + arv1.divideByLeft (arv2));
		// System.out.println ("dotMult: " + arv1.dotMult (arv2));
	}
}
// end of file
