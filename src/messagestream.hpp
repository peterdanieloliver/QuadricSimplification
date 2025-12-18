#pragma once

#include <qtextedit>
#include <iostream>
#include <fstream>
#include <string>

// class for redirecting std::cout buffer to custom buffer that can be read by the textEdit widget
class MessageStream : public std::basic_streambuf<char>
{

private:

	std::ostream& m_stream;
	std::streambuf* m_old_buf;
	std::string m_string;
	QTextEdit* message_box;

public:

	MessageStream(std::ostream& stream, QTextEdit* text_edit)
		:m_stream(stream)
	{
		this->message_box = text_edit;
		this->m_old_buf = stream.rdbuf();

		stream.rdbuf(this);
	}

	virtual ~MessageStream()
	{
		if (!m_string.empty())
		{
			message_box->append(m_string.c_str());
		}

		this->m_stream.rdbuf(this->m_old_buf);
	}

protected:

	virtual int_type overflow(int_type v)
	{
		if (v == '\n')
		{
			message_box->append(m_string.c_str());
			m_string.erase(m_string.begin(), m_string.end());
		}
		else
		{
			return v;
		}
	}

	virtual std::streamsize xsputn(const char* p, std::streamsize n)
	{
		m_string.append(p, p + n);
		int pos = 0;
		while (pos != std::string::npos)
		{
			pos = m_string.find('\n');
			if (pos != std::string::npos)
			{
				std::string tmp(m_string.begin(), m_string.end() + pos);
				message_box->append(tmp.c_str());
				m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
			}
		}
		return n;
	}

};

void printLine(std::string line, bool new_section = false)
{
	if (new_section) { std::cout << std::endl; }
	std::cout << line.c_str() << std::endl;
}
